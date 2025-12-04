import torch
from tqdm import tqdm
import torch.nn.functional as F
from .model import MEDR

def target_distribution(batch):
    weight = (batch ** 2) / torch.sum(batch, 0)
    return (weight.t() / torch.sum(weight, 1)).t()


def reconstruction_loss(decoded, x):
    loss_func = torch.nn.MSELoss()
    loss_rcn = loss_func(decoded, x)
    return loss_rcn


# def gcn_loss(preds, labels, mu, logvar, n_nodes, norm, mask=None):
#     if mask is not None:
#         preds = preds * mask
#         labels = labels * mask
#
#     cost = norm * F.binary_cross_entropy_with_logits(preds, labels)
#
#     # see Appendix B from VAE paper:
#     # Kingma and Welling. Auto-Encoding Variational Bayes. ICLR, 2014
#     # https://arxiv.org/abs/1312.6114
#     # 0.5 * sum(1 + log(sigma^2) - mu^2 - sigma^2)
#     KLD = -0.5 / n_nodes * torch.mean(torch.sum(
#         1 + 2 * logvar - mu.pow(2) - logvar.exp().pow(2), 1))
#     return cost + KLD



def gcn_loss(preds, labels, mu, logvar, n_nodes, norm):
    cost = norm * F.binary_cross_entropy_with_logits(preds, labels)
    # see Appendix B from VAE paper:
    # Kingma and Welling. Auto-Encoding Variational Bayes. ICLR, 2014
    # https://arxiv.org/abs/1312.6114
    # 0.5 * sum(1 + log(sigma^2) - mu^2 - sigma^2)
    KLD = -0.5 / n_nodes * torch.mean(torch.sum(
        1 + 2 * logvar - mu.pow(2) - logvar.exp().pow(2), 1))
    return cost + KLD


class Trainer:
    def __init__(
            self,
            X,
            graph_spatial,
            graph_feat,
            device = 'cuda:0',
            dim_input = 3000,
            dim_output = 32
    ):

        self.cell_num = len(X)

        self.device = device if torch.cuda.is_available() else 'cpu'

        print(f"Using device: {self.device}")

        self.X = torch.FloatTensor(X.copy()).to(self.device)
        self.dim_input = self.X.shape[1]
        self.dim_output = dim_output

        self.adj_spatial_norm = graph_spatial["adj_norm"].to(self.device)
        self.adj_spatial_label = graph_spatial["adj_label"].to(self.device)
        self.adj_spatial_norm_value = graph_spatial["norm_value"]

        self.adj_feat_norm = graph_feat["adj_norm"].to(self.device)
        self.adj_feat_label = graph_feat["adj_label"].to(self.device)
        self.adj_feat_norm_value = graph_feat["norm_value"]

        # self.model = MEDR(dim_input=self.dim_input, dim_output=self.dim_output).to(self.device)
        # self.model = MEDR(in_dim=self.dim_input, num_hidden=128, out_dim=32).to(self.device)
        self.model = MEDR(self.dim_input).to(self.device)

        self.mask = False


    def mask_generator(self, adj_label, N=1):
        idx = adj_label.indices()

        list_non_neighbor = []
        for i in range(0, self.cell_num):
            neighbor = idx[1, torch.where(idx[0, :] == i)[0]]
            n_selected = len(neighbor) * N

            # non neighbors
            total_idx = torch.range(0, self.cell_num-1, dtype=torch.float32).to(self.device)
            non_neighbor = total_idx[~torch.isin(total_idx, neighbor)]
            indices = torch.randperm(len(non_neighbor), dtype=torch.float32).to(self.device)
            random_non_neighbor = indices[:n_selected]
            list_non_neighbor.append(random_non_neighbor)

        x = torch.repeat_interleave(adj_label.indices()[0], N)
        y = torch.concat(list_non_neighbor)

        indices = torch.stack([x, y])
        indices = torch.concat([adj_label.indices(), indices], axis=1)

        value = torch.concat([adj_label.values(), torch.zeros(len(x), dtype=torch.float32).to(self.device)])
        adj_mask = torch.sparse_coo_tensor(indices, value)

        return adj_mask
    

    def train(
            self,
            epochs=200,
            lr=0.01,
            decay=0.01
    ):
        self.optimizer = torch.optim.Adam(
            params=list(self.model.parameters()),
            lr=lr,
            weight_decay=decay)

        self.model.train()

        for epoch in tqdm(range(epochs)):
            self.model.train()
            self.optimizer.zero_grad()


            z, recon, pre_recon, loss_self, pre_loss_self, mu, logvar, pre_mu, pre_logvar, z1 = self.model(self.X, self.adj_spatial_norm, self.adj_feat_norm)

            if self.mask:
                pass
            else:

                adj_mask = self.mask_generator(self.adj_spatial_label, N=1)
                self.adj_mask = adj_mask
                
                # pre_adj_mask = self.mask_generator(self.adj_feat_label, N=1)
                # self.pre_adj_mask = pre_adj_mask

                self.mask = True
            
            loss_gcn = gcn_loss(
                preds=self.model.dc(z, self.adj_mask),
                # labels=self.adj_label,
                labels=self.adj_mask.coalesce().values(),
                mu=mu,
                logvar=logvar,
                n_nodes=self.cell_num,
                norm=self.adj_spatial_norm_value
            )

            # pre_loss_gcn = gcn_loss(
            #     preds=self.model.dc(z, self.pre_adj_mask),
            #     # labels=self.adj_label,
            #     labels=self.pre_adj_mask.coalesce().values(),
            #     mu=pre_mu,
            #     logvar=pre_logvar,
            #     n_nodes=self.cell_num,
            #     norm=self.adj_feat_norm_value
            # )




            loss_rec = reconstruction_loss(recon, self.X)

            pre_loss_rec = reconstruction_loss(pre_recon, self.X)
           
            

            # loss = 10*reconstruction_loss(recon, self.X) + 0.5*reconstruction_loss(emb_latent_recon, emb_latent_feature)
            # loss = 10 * reconstruction_loss(recon, self.X) + 1 * loss_self

            # loss = 10*loss_rec + pre_loss_rec + 0.1*loss_gcn + 0.1*pre_loss_gcn + loss_self + pre_loss_self
            # loss = 10*loss_rec + loss_gcn + pre_loss_gcn + 5*loss_self + 5*pre_loss_rec
            #10 5 1
            loss = 7*loss_rec + 3*loss_gcn + 1*loss_self
            # loss = 10*loss_rec + 5*loss_gcn + 1*loss_self

            # loss = 30*pre_loss_rec
            loss.backward() 
            self.optimizer.step() 


    def process(self):
        self.model.eval()
        z, recon, pre_recon, loss_self, pre_loss_self, mu, logvar, pre_mu, pre_logvar, z1 = self.model(self.X, self.adj_spatial_norm, self.adj_feat_norm)
        
        latent = z.data.cpu().numpy()

        latent1 = z1.data.cpu().numpy()

        
        return latent, latent1