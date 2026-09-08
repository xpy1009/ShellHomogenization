import torch
import numpy as np
import matplotlib.pyplot as plt
import time

class Model(torch.nn.Module):
    def __init__(self, scale, input_dimension, output_dimension, n_hidden_layers, neurons):
        super(Model, self).__init__()
        self.scale = 1.0 / scale
        self.activation1 = torch.nn.SiLU()
        self.activation2 = torch.nn.Softplus()
        self.input_layer = torch.nn.Linear(input_dimension, neurons)
        self.hidden_layers = torch.nn.ModuleList([torch.nn.Linear(neurons, neurons) for _ in range(n_hidden_layers - 1)])
        self.output_layer = torch.nn.Linear(neurons, output_dimension)

    def forward(self, x):
        x = self.activation1(self.input_layer(x))
        for _, l in enumerate(self.hidden_layers):
            x = self.activation1(l(x))
        x = self.activation2(self.output_layer(x))
        return x
    
    @torch.jit.export
    def energy(self, x: torch.Tensor) -> torch.Tensor:
        x.requires_grad_(True)
        y: torch.Tensor = self.scale * self.forward(x)
        assert y is not None
        return y
    
    @torch.jit.export
    def gradient(self, x: torch.Tensor) -> torch.Tensor:
        x.requires_grad_(True)
        y: torch.Tensor = self.energy(x)
        grads_y: Optional[torch.Tensor] = torch.ones_like(y)
        dydx = torch.autograd.grad([y], [x], [grads_y], create_graph=True)[0]
        assert dydx is not None
        return dydx
    
    @torch.jit.export
    def hessian(self, x: torch.Tensor):
        dydx = self.gradient(x)
        dydx = dydx.sum(0)
        grads: Optional[torch.Tensor] = torch.ones_like(dydx[0])
        H_row0 = torch.autograd.grad([dydx[0]], [x], [grads], create_graph=True)[0]
        H_row1 = torch.autograd.grad([dydx[1]], [x], [grads], create_graph=True)[0]
        H_row2 = torch.autograd.grad([dydx[2]], [x], [grads], create_graph=True)[0]
        assert H_row0 is not None
        assert H_row1 is not None
        assert H_row2 is not None
        return torch.stack((H_row0, H_row1, H_row2),1)

def energy_loss(pred, gt):
    pred_normalized = pred / gt
    gt_normalized = gt / gt
    loss = gt_normalized - pred_normalized
    return loss

def moment_loss(pred, gt):
    moment_norm = torch.norm(gt, dim=1, keepdim=True)
    pred_normalized = pred / moment_norm
    gt_normalized = gt / moment_norm
    loss = gt_normalized - pred_normalized
    return loss
    
def o3_loss(model, n_samples):
    x1 = 400 * (torch.rand(n_samples, dtype=torch.float64, requires_grad=True) - 0.5)
    x2 = 400 * (torch.rand(n_samples, dtype=torch.float64, requires_grad=True) - 0.5)
    t = np.pi * torch.rand(n_samples, dtype=torch.float64)
    x = torch.stack([x1*torch.cos(t)**2+x2*torch.sin(t)**2,
                x1*torch.sin(t)**2+x2*torch.cos(t)**2,
                2*(x1-x2)*torch.cos(t)*torch.sin(t)], 1)

    y = model(x)

    dydx = torch.autograd.grad(y.sum(), (x1,x2), create_graph=True)
    dydx1 = dydx[0]
    dydx2 = dydx[1]

    d2ydx12 = torch.autograd.grad(dydx1.sum(), x1, create_graph=True)[0]
    d2ydx22 = torch.autograd.grad(dydx2.sum(), x2, create_graph=True)[0]

    d3ydx12dxi = torch.autograd.grad(d2ydx12.sum(), (x1,x2), create_graph=True)
    d3ydx13 = d3ydx12dxi[0]
    d3ydx12dx2 = d3ydx12dxi[1]
    d3ydx22dxi = torch.autograd.grad(d2ydx22.sum(), (x1,x2), create_graph=True)
    d3ydx1dx22 = d3ydx22dxi[0]
    d3ydx23 = d3ydx22dxi[1]

    loss = d3ydx13.square()+d3ydx12dx2.square()+d3ydx1dx22.square()+d3ydx23.square()
    return loss


def fit(training_loader, model, num_epochs, optimizer, aqap, verbose):
    history = list()
    o3_weight = 1e2

    # Loop over epochs
    for epoch in range(num_epochs):
        epoch_history = list()
        for _, (kj_, mj_, psi_) in enumerate(training_loader):
            def closure():
                kj, mj, psi = kj_, mj_, psi_
                optimizer.zero_grad()
                kj.requires_grad = True
                psi_pred = model(kj)
                loss = torch.mean(torch.square(energy_loss(psi_pred, psi)))

                m_pred = torch.autograd.grad(psi_pred.sum(), kj, create_graph=True)[0]
                loss += torch.mean(torch.square(moment_loss(m_pred, mj)))

                if (aqap):
                    loss += o3_weight * torch.mean(o3_loss(model, 64))

                loss.backward()

                epoch_history.append(loss.item())
                if (torch.isnan(loss).any()):
                    assert(0)
                    
                return loss

            optimizer.step(closure=closure)
        losses = 0.
        for loss in epoch_history:
            losses += loss 
        losses = losses / len(epoch_history)
        history.append(losses)

        if (verbose):
            print(epoch, ": ", losses)

    if (verbose):
        print('Final Loss: ', history[-1])

    return history


def train(data, scale, name=None, aqap=False, verbose=False):
    es = data[:,:3]
    ss = data[:,3:6]
    psis = data[:,6:7]

    idx = (psis > 1e-5).squeeze(1)
    es = es[idx,:]
    ss = ss[idx,:]
    psis = psis[idx]

    if (verbose):
        # check data scale
        print(psis.max(), ss.max())
        print(psis.min(), ss.norm(dim=1).abs().min())
        print(es.shape)

    # prepare data
    generator = torch.Generator().manual_seed(42)
    training_set, validation_set = torch.utils.data.random_split(torch.utils.data.TensorDataset(es, scale * ss, scale * psis), [.95, .05], generator=generator)

    batch_size = 16
    training_loader = torch.utils.data.DataLoader(training_set, batch_size=batch_size, shuffle=True)
    validation_loader = torch.utils.data.DataLoader(validation_set, batch_size=batch_size, shuffle=False)

    # model
    model = Model(scale, 3, 1, 4, 32)
    model.double()

    # training
    n_epochs = 1000
    optimizer_ADAM = torch.optim.Adam(model.parameters(), lr=float(0.001))
    history = fit(training_loader, model, n_epochs, optimizer_ADAM, aqap, verbose)

    n_epochs = 500
    optimizer_ADAM = torch.optim.Adam(model.parameters(), lr=float(0.0001))
    history += fit(training_loader, model, n_epochs, optimizer_ADAM, aqap, verbose)

    # log
    if (name is not None):
        fig, ax = plt.subplots()
        ax.set_yscale('log')
        ax.plot(history)
        fig.savefig(name + '.png')

        es.requires_grad = True
        psi_pred = model.energy(es)
        s_pred = model.gradient(es)
        psi_err = energy_loss(psi_pred, psis)
        s_err = moment_loss(s_pred, ss)
        es.requires_grad = False
        if (verbose):
            print(psi_err.max().item(), s_err.max().item(), o3_loss(model, es.shape[0]).max().item())

    return model



def main():
    start = time.time()

    verbose = False

    print("train stretching model...")
    scale = 1e0
    torch.manual_seed(42)
    stretch_data = torch.from_numpy(np.loadtxt("../data/stretching.txt"))[:, 3:]
    stretch_model = train(stretch_data, scale, "../data/stretch", False, verbose)
    torch.jit.script(stretch_model).save('../data/stretch.pt')

    print("train bending model...")
    torch.manual_seed(42)
    scale = 1e2
    bend_data = torch.from_numpy(np.loadtxt("../data/bending.txt"))[:, 2:]
    bend_model = train(bend_data, scale, "../data/bend", True, verbose)
    torch.jit.script(bend_model).save('../data/bend.pt')

    end = time.time()
    print(end - start, "s")


if __name__ == "__main__":
    main()
