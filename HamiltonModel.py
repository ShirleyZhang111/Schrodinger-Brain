import torch
import torch.nn as nn
import os
os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'  

# Hamilton Model: learns auxiliary signals and linear Hamiltonian H from observed signals
class Hamilton_Model(nn.Module):
    def __init__(self, m):
        """
        Initialize Hamilton Model
        
        Args:
            m (int): Dimension of the system (number of regions)
        """
        super(Hamilton_Model, self).__init__()
        self.m = m

        # Hamiltonian matrix parameters
        self.log_diag = nn.Parameter(torch.randn(m))      # Logarithm of diagonal elements (ensures positivity)
        self.off_diag = nn.Parameter(torch.randn(m, m))   # Off-diagonal elements

         # CNN for learning auxiliary signals y from observed signals x
        self.f_net = nn.Sequential(
            # Input: (batch_size=m, 1, num_samples)
            nn.Conv1d(in_channels=1, out_channels=16, kernel_size=3, padding=1),
            nn.ReLU(),
            nn.Conv1d(in_channels=16, out_channels=32, kernel_size=3, padding=1),
            nn.ReLU(),
            nn.Conv1d(in_channels=32, out_channels=16, kernel_size=3, padding=1),
            nn.ReLU(),
            nn.Conv1d(in_channels=16, out_channels=1, kernel_size=3, padding=1),
        )
    
    def forward(self, x):
        """
        Forward pass of the Hamilton Model
        
        Args:
            x (torch.Tensor): Input observed signals of shape (num_samples, m)
            
        Returns:
            H (torch.Tensor): Hamiltonian matrix of shape (m, m)
            y (torch.Tensor): Auxiliary signals of shape (num_samples, m)
        """

        # Construct symmetric Hamiltonian matrix H
        D = torch.diag(torch.exp(self.log_diag))      # Diagonal matrix with positive elements
        L = torch.tril(self.off_diag, diagonal=-1)    # Lower triangular matrix (excluding diagonal)
        H = D + 0.5*(L + L.T)                         # Symmetric matrix: H = D + (L + L^T)/2
        
       
        # (num_samples, m) -> (m, 1, num_samples)
        x_cnn = x.permute(1, 0).unsqueeze(1)  
               
        # y: (1, num_samples, m)
        y_cnn = self.f_net(x_cnn)

        # (m, 1, num_samples) -> (num_samples, m)
        y = y_cnn.squeeze(1).permute(1, 0)        

        return H, y
    
    def compute_loss(self, x, dx_dt_obs,dt):
        """
        Compute loss function based on Hamilton's equations
        
        Args:
            x (torch.Tensor): Observed signals of shape (num_samples, m)
            dx_dt_obs (torch.Tensor): Observed derivatives of x of shape (num_samples-1, m)
            dt (float): Time step size
            
        Returns:
            total_loss (torch.Tensor): Combined loss from both Hamilton equations
            loss1 (torch.Tensor): Loss from first Hamilton equation
            loss2 (torch.Tensor): Loss from second Hamilton equation
        """
        H, y = self.forward(x)  # y: (num_samples, m)
        
        # Hamilton's equation
        # dx/dt = -Hy
        dx_dt_pred = -torch.einsum("ij,sj->si", H, y[:-1])  # -H @ y[:-1].T
        loss1 = torch.mean((dx_dt_pred - dx_dt_obs) ** 2)
        
        # dy/dt = Hx
        dy_dt_pred = (y[1:] - y[:-1]) / dt  
        dy_dt_obs = torch.einsum("ij,sj->si", H, x[:-1])  # H @ x[:-1].T
        loss2 = torch.mean((dy_dt_pred - dy_dt_obs) ** 2)
        
        total_loss = loss1 + loss2
        return total_loss, loss1, loss2