import numpy as np
import torch
from scipy import stats
import matplotlib.pyplot as plt
from scipy.signal import hilbert
import os
from scipy.io import loadmat
from HamiltonModel import Hamilton_Model
import torch.optim as optim
os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'

# Load regional resting-state fMRI data 
dict_mat = loadmat(".\\Data\\TC_HCP_379.mat")
x = dict_mat["TC"]  # shape: (spatial, time)
x = stats.zscore(x, axis=1) # z-scored
x = x[:,100:1000] # use 101~1000 time points (length = 900 TRs)

m = x.shape[1] # Number of spatial regions
num_samples = x.shape[0] # Number of time points (TRs)

x = torch.from_numpy(x).float() 
dt = 1 # Time step (TR)


# Compute empirical derivative using finite differences
dx_dt_obs = torch.zeros(num_samples-1, m)
for t in range(num_samples-1):
    dx_dt_obs[t] = (x[t+1] - x[t]) / dt
dx_dt_obs = dx_dt_obs.float()


# Initialize Hamilton neural network model and optimizer
model = Hamilton_Model(m)
optimizer = optim.Adam(model.parameters(), lr=0.001)

# Model training configuration
num_epochs = 10000
loss_history = []

# Training loop
for epoch in range(num_epochs):
    optimizer.zero_grad()
    total_loss, loss1, loss2 = model.compute_loss(x, dx_dt_obs,dt)
    total_loss.backward()
    optimizer.step()   
    loss_history.append(total_loss.item())
    
    if epoch % 100 == 0:
        print(f'Epoch {epoch}, Loss: {total_loss.item():.4f}, Loss1: {loss1.item():.4f}, Loss2: {loss2.item():.4f}')


# Save trained model parameters
torch.save(model, '.\\Results\\HamiltonModel_Conv1d_900TRs1.pt')
model.eval()


# Apply discrete Hilbert transform to rs-fMRI signal 'x'
x = dict_mat["Data"]
x_new = np.matrix.transpose(x)
hilbert_transformed = hilbert(x_new, axis=0)
hilbert_imag = np.imag(hilbert_transformed)  
Hilbert_y = torch.from_numpy(hilbert_imag).float()
x_new = torch.from_numpy(x_new).float()

# Comparison between Learned Auxiliary Signal and Discrete Hilbert Transform
with torch.no_grad():
    H_learned, y_learned = model(x_new)
    plt.figure(figsize=(10, 4))

    # Plot examples for 4 spatial regions
    for i in range(4): 
        plt.subplot(2, 2, i+1)
        plt.plot(Hilbert_y[101:300, i], label = 'Discrete Hilbert Transform' )
        plt.plot(y_learned[101:300, i], label='Learned Auxiliary Signal')
        plt.title(f'Dimension {i+1}')

        # Calculate and display correlation coefficient
        correlation = np.corrcoef(Hilbert_y[101:300, i], y_learned[101:300, i])[0, 1]
        plt.text(0.05, 0.95, f'Correlation: {correlation:.3f}', 
                transform=plt.gca().transAxes,
                ha='left', va='top',
                bbox=dict(facecolor='white', alpha=0.8))         
        if i == 1:
            plt.legend(loc='best')

    plt.suptitle('Comparison: Hilbert Transform vs Learned Auxiliary Signal')
    plt.show()


# Save Learned Parameters (Auxilaray signal y and coupling matrix H)
H_learned = H_learned.detach().numpy()
y_learned = y_learned.detach().numpy()
np.save('.\\Results\\H_learned.npy', H_learned)
np.save('.\\Results\\y_learned.npy', y_learned)