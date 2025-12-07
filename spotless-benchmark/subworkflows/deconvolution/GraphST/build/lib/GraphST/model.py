import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.nn.parameter import Parameter
from torch.nn.modules.module import Module
from torch_geometric.nn import GATConv, SGConv, GATv2Conv

class Discriminator(nn.Module):
    def __init__(self, n_h):
        super(Discriminator, self).__init__()
        self.f_k = nn.Bilinear(n_h, n_h, 1)

        for m in self.modules():
            self.weights_init(m)

    def weights_init(self, m):
        if isinstance(m, nn.Bilinear):
            torch.nn.init.xavier_uniform_(m.weight.data)
            if m.bias is not None:
                m.bias.data.fill_(0.0)

    def forward(self, c, h_pl, h_mi, s_bias1=None, s_bias2=None):
        c_x = c.expand_as(h_pl)  

        sc_1 = self.f_k(h_pl, c_x)
        sc_2 = self.f_k(h_mi, c_x)

        if s_bias1 is not None:
            sc_1 += s_bias1
        if s_bias2 is not None:
            sc_2 += s_bias2

        logits = torch.cat((sc_1, sc_2), 1)

        return logits
    
class AvgReadout(nn.Module):
    def __init__(self):
        super(AvgReadout, self).__init__()

    def forward(self, emb, mask=None):
        vsum = torch.mm(mask, emb)
        row_sum = torch.sum(mask, 1)
        row_sum = row_sum.expand((vsum.shape[1], row_sum.shape[0])).T
        global_emb = vsum / row_sum 
          
        return F.normalize(global_emb, p=2, dim=1) 
    
class Encoder(Module):
    def __init__(self, in_features, out_features, graph_neigh, dropout=0.0, act=F.relu):
        super(Encoder, self).__init__()
        self.in_features = in_features
        self.out_features = out_features
        self.graph_neigh = graph_neigh
        self.dropout = dropout
        self.act = act
        
        self.weight1 = Parameter(torch.FloatTensor(self.in_features, self.out_features))
        self.weight2 = Parameter(torch.FloatTensor(self.out_features, self.in_features))
        self.reset_parameters()
        
        self.disc = Discriminator(self.out_features)

        self.sigm = nn.Sigmoid()
        self.read = AvgReadout()
        
    def reset_parameters(self):
        torch.nn.init.xavier_uniform_(self.weight1)
        torch.nn.init.xavier_uniform_(self.weight2)

    def forward(self, feat, feat_a, adj):
        z = F.dropout(feat, self.dropout, self.training)
        z = torch.mm(z, self.weight1)
        z = torch.mm(adj, z)
        
        hiden_emb = z
        
        h = torch.mm(z, self.weight2)
        h = torch.mm(adj, h)
        
        emb = self.act(z)
        
        z_a = F.dropout(feat_a, self.dropout, self.training)
        z_a = torch.mm(z_a, self.weight1)
        z_a = torch.mm(adj, z_a)
        emb_a = self.act(z_a)
        
        g = self.read(emb, self.graph_neigh) 
        g = self.sigm(g)  

        g_a = self.read(emb_a, self.graph_neigh)
        g_a = self.sigm(g_a)  

        ret = self.disc(g, emb, emb_a)  
        ret_a = self.disc(g_a, emb_a, emb) 
        
        return hiden_emb, h, ret, ret_a

class Encoder_GAT(Module):
    def __init__(self, in_features, out_features, graph_neigh, dropout=0.0, act=F.relu):
        super(Encoder_GAT, self).__init__()
        self.in_features = in_features
        self.out_features = out_features
        self.graph_neigh = graph_neigh
        self.dropout = dropout
        self.act = act
        
        # Replace manual weights with GAT layers
        # Heads=1 keeps dimensions simple (in -> out). 
        # You can increase heads but must adjust concat/out_channels accordingly.
        self.conv1 = GATConv(self.in_features, self.out_features, heads=1, concat=True, dropout=dropout)
        self.conv2 = GATConv(self.out_features, self.in_features, heads=1, concat=True, dropout=dropout)

        self.disc = Discriminator(self.out_features)
        self.sigm = nn.Sigmoid()
        self.read = AvgReadout()

    def forward(self, feat, feat_a, adj):
        # Convert the Adjacency Matrix (adj) to Edge Index (COO format)
        # GraphST passes 'adj' as a dense tensor or sparse tensor, but GAT needs edge indices.
        if adj.is_sparse:
            edge_index = adj._indices()
        else:
            edge_index = adj.nonzero().t().contiguous()

        # --- Branch 1: Main Features ---
        # Original: z = torch.mm(z, self.weight1); z = torch.mm(adj, z)
        # GAT: Input features -> Attention-weighted aggregation
        z = self.conv1(feat, edge_index)
        z = F.dropout(z, self.dropout, self.training)
        
        hiden_emb = z
        
        # --- Reconstruction (Decoder step) ---
        # Original: h = torch.mm(z, self.weight2); h = torch.mm(adj, h)
        h = self.conv2(z, edge_index)
        
        emb = self.act(z)
        
        # --- Branch 2: Augmented Features (Contrastive View) ---
        # Shared weights (conv1) used for augmented view
        z_a = self.conv1(feat_a, edge_index)
        z_a = F.dropout(z_a, self.dropout, self.training)
        emb_a = self.act(z_a)
        
        # --- Readout & Discriminator (Unchanged) ---
        g = self.read(emb, self.graph_neigh) 
        g = self.sigm(g)  

        g_a = self.read(emb_a, self.graph_neigh)
        g_a = self.sigm(g_a)  

        ret = self.disc(g, emb, emb_a)  
        ret_a = self.disc(g_a, emb_a, emb) 
        
        return hiden_emb, h, ret, ret_a

class Encoder_GATv2(Module):
    def __init__(self, in_features, out_features, graph_neigh, dropout=0.0, act=F.relu):
        super(Encoder_GATv2, self).__init__()
        self.in_features = in_features
        self.out_features = out_features
        self.graph_neigh = graph_neigh
        self.dropout = dropout
        self.act = act
        
        # GATv2Conv: Fixes static attention of standard GAT
        # heads=1, concat=True keeps output dim = out_features (if heads>1, output is heads*out)
        self.conv1 = GATv2Conv(self.in_features, self.out_features, heads=1, concat=True, dropout=dropout)
        self.conv2 = GATv2Conv(self.out_features, self.in_features, heads=1, concat=True, dropout=dropout)

        self.disc = Discriminator(self.out_features)
        self.sigm = nn.Sigmoid()
        self.read = AvgReadout()

    def forward(self, feat, feat_a, adj):
        # Convert adj to edge_index (required for PyG layers)
        if adj.is_sparse:
            edge_index = adj._indices()
        else:
            edge_index = adj.nonzero().t().contiguous()

        # --- Branch 1: Main Features ---
        z = self.conv1(feat, edge_index)
        z = F.dropout(z, self.dropout, self.training)
        
        hiden_emb = z
        
        # Decoder
        h = self.conv2(z, edge_index)
        
        emb = self.act(z)
        
        # --- Branch 2: Augmented Features ---
        z_a = self.conv1(feat_a, edge_index)
        z_a = F.dropout(z_a, self.dropout, self.training)
        emb_a = self.act(z_a)
        
        # --- Readout & Discriminator ---
        g = self.read(emb, self.graph_neigh) 
        g = self.sigm(g)  

        g_a = self.read(emb_a, self.graph_neigh)
        g_a = self.sigm(g_a)  

        ret = self.disc(g, emb, emb_a)  
        ret_a = self.disc(g_a, emb_a, emb) 
        
        return hiden_emb, h, ret, ret_a


class Encoder_SGC(Module):
    def __init__(self, in_features, out_features, graph_neigh, dropout=0.0, act=F.relu):
        super(Encoder_SGC, self).__init__()
        self.in_features = in_features
        self.out_features = out_features
        self.graph_neigh = graph_neigh
        self.dropout = dropout
        
        # SGC Layer: K=2 aggregates information from 2-hop neighbors in one step
        # cached=False ensures it works with changing graphs (if any), though usually True is faster for static graphs.
        self.conv1 = SGConv(self.in_features, self.out_features, K=2, cached=False)
        self.conv2 = SGConv(self.out_features, self.in_features, K=2, cached=False)

        self.disc = Discriminator(self.out_features)
        self.sigm = nn.Sigmoid()
        self.read = AvgReadout()

    def forward(self, feat, feat_a, adj):
        # Convert adj to edge_index (required for PyG layers)
        if adj.is_sparse:
            edge_index = adj._indices()
        else:
            edge_index = adj.nonzero().t().contiguous()

        # --- Branch 1: Main Features ---
        z = F.dropout(feat, self.dropout, self.training)
        z = self.conv1(z, edge_index)
        
        hiden_emb = z
        
        # Reconstruction (Decoder)
        h = self.conv2(z, edge_index)
        
        # SGC typically removes non-linearities, so we skip self.act(z)
        emb = z 
        
        # --- Branch 2: Augmented Features ---
        z_a = F.dropout(feat_a, self.dropout, self.training)
        z_a = self.conv1(z_a, edge_index)
        emb_a = z_a # Skip activation
        
        # --- Readout & Discriminator ---
        g = self.read(emb, self.graph_neigh) 
        g = self.sigm(g)  

        g_a = self.read(emb_a, self.graph_neigh)
        g_a = self.sigm(g_a)  

        ret = self.disc(g, emb, emb_a)  
        ret_a = self.disc(g_a, emb_a, emb) 
        
        return hiden_emb, h, ret, ret_a

class Encoder_sparse(Module):
    """
    Sparse version of Encoder
    """
    def __init__(self, in_features, out_features, graph_neigh, dropout=0.0, act=F.relu):
        super(Encoder_sparse, self).__init__()
        self.in_features = in_features
        self.out_features = out_features
        self.graph_neigh = graph_neigh
        self.dropout = dropout
        self.act = act
        
        self.weight1 = Parameter(torch.FloatTensor(self.in_features, self.out_features))
        self.weight2 = Parameter(torch.FloatTensor(self.out_features, self.in_features))
        self.reset_parameters()
        
        self.disc = Discriminator(self.out_features)

        self.sigm = nn.Sigmoid()
        self.read = AvgReadout()
        
    def reset_parameters(self):
        torch.nn.init.xavier_uniform_(self.weight1)
        torch.nn.init.xavier_uniform_(self.weight2)

    def forward(self, feat, feat_a, adj):
        z = F.dropout(feat, self.dropout, self.training)
        z = torch.mm(z, self.weight1)
        z = torch.spmm(adj, z)
        
        hiden_emb = z
        
        h = torch.mm(z, self.weight2)
        h = torch.spmm(adj, h)
        
        emb = self.act(z)
        
        z_a = F.dropout(feat_a, self.dropout, self.training)
        z_a = torch.mm(z_a, self.weight1)
        z_a = torch.spmm(adj, z_a)
        emb_a = self.act(z_a)
         
        g = self.read(emb, self.graph_neigh)
        g = self.sigm(g)
        
        g_a = self.read(emb_a, self.graph_neigh)
        g_a =self.sigm(g_a)       
       
        ret = self.disc(g, emb, emb_a)  
        ret_a = self.disc(g_a, emb_a, emb)
        
        return hiden_emb, h, ret, ret_a     

class Encoder_sc(torch.nn.Module):
    def __init__(self, dim_input, dim_output, dropout=0.0, act=F.relu):
        super(Encoder_sc, self).__init__()
        self.dim_input = dim_input
        self.dim1 = 256
        self.dim2 = 64
        self.dim3 = 32
        self.act = act
        self.dropout = dropout
        
        #self.linear1 = torch.nn.Linear(self.dim_input, self.dim_output)
        #self.linear2 = torch.nn.Linear(self.dim_output, self.dim_input)
        
        #self.weight1_en = Parameter(torch.FloatTensor(self.dim_input, self.dim_output))
        #self.weight1_de = Parameter(torch.FloatTensor(self.dim_output, self.dim_input))
        
        self.weight1_en = Parameter(torch.FloatTensor(self.dim_input, self.dim1))
        self.weight2_en = Parameter(torch.FloatTensor(self.dim1, self.dim2))
        self.weight3_en = Parameter(torch.FloatTensor(self.dim2, self.dim3))
        
        self.weight1_de = Parameter(torch.FloatTensor(self.dim3, self.dim2))
        self.weight2_de = Parameter(torch.FloatTensor(self.dim2, self.dim1))
        self.weight3_de = Parameter(torch.FloatTensor(self.dim1, self.dim_input))
      
        self.reset_parameters()

    def reset_parameters(self):
        torch.nn.init.xavier_uniform_(self.weight1_en)
        torch.nn.init.xavier_uniform_(self.weight1_de)
        
        torch.nn.init.xavier_uniform_(self.weight2_en)
        torch.nn.init.xavier_uniform_(self.weight2_de)
        
        torch.nn.init.xavier_uniform_(self.weight3_en)
        torch.nn.init.xavier_uniform_(self.weight3_de)
        
    def forward(self, x):
        x = F.dropout(x, self.dropout, self.training)
        
        #x = self.linear1(x)
        #x = self.linear2(x)
        
        #x = torch.mm(x, self.weight1_en)
        #x = torch.mm(x, self.weight1_de)
        
        x = torch.mm(x, self.weight1_en)
        x = torch.mm(x, self.weight2_en)
        x = torch.mm(x, self.weight3_en)
        
        x = torch.mm(x, self.weight1_de)
        x = torch.mm(x, self.weight2_de)
        x = torch.mm(x, self.weight3_de)
        
        return x
    
class Encoder_map(torch.nn.Module):
    def __init__(self, n_cell, n_spot):
        super(Encoder_map, self).__init__()
        self.n_cell = n_cell
        self.n_spot = n_spot
          
        self.M = Parameter(torch.FloatTensor(self.n_cell, self.n_spot))
        self.reset_parameters()

    def reset_parameters(self):
        torch.nn.init.xavier_uniform_(self.M)
        
    def forward(self):
        x = self.M
        
        return x 
    

