import numpy as np
from sklearn.decomposition import PCA
from rdkit import Chem
from rdkit.Chem import Descriptors
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd

# Define a function to calculate molecular descriptors
def calc_descriptors(mol):
    return [Descriptors.MolWt(mol), Descriptors.TPSA(mol), Descriptors.NumRotatableBonds(mol)]

# Load the data
location = '/Users/ronannovosel/Downloads/logP_dataset.csv'
df = pd.read_csv(location, header=None)
smiles = df[0].tolist()
logP = df[1].tolist()

# Convert SMILES strings to molecules and check for invalid SMILES
mols = []
valid_logP = []
invalid_smiles = []

for i, smi in enumerate(smiles):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        invalid_smiles.append(smi)
    else:
        mols.append(mol)
        valid_logP.append(logP[i])

# Print invalid SMILES
if invalid_smiles:
    print('Invalid SMILES strings:')
    for smi in invalid_smiles:
        print(smi)

# Calculate molecular descriptors using list comprehension
X = np.array([calc_descriptors(mol) for mol in mols])

# Perform PCA to reduce the dimensionality of the data
pca = PCA(n_components=3)
X_pca = pca.fit_transform(X)

# Create a 3D plot with black background and adjusted color map
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection='3d')
scatter = ax.scatter(X_pca[:, 0], X_pca[:, 1], X_pca[:, 2], c=valid_logP, cmap='plasma', alpha=0.7)
ax.set_facecolor('black')
ax.set_xlabel('PC1', fontsize=12, labelpad=12)
ax.set_ylabel('PC2', fontsize=12, labelpad=12)
ax.set_zlabel('PC3', fontsize=12, labelpad=12)
ax.tick_params(axis='x', labelsize=10)
ax.tick_params(axis='y', labelsize=10)
ax.tick_params(axis='z', labelsize=10)

# Customize colorbar
cbar = plt.colorbar(scatter, fraction=0.046, pad=0.04)
cbar.set_label('LogP', fontsize=12)
cbar.ax.tick_params(labelsize=10)

# Add a title and adjust plot layout
plt.title('Chemical Space', fontsize=14, pad=20)
plt.tight_layout()

# Add a legend to show which molecular descriptors were used
legend_elements = [
    plt.Line2D([], [], marker='o', color='w', markerfacecolor='black', markersize=8, label='Molecular Weight'),
    plt.Line2D([], [], marker='o', color='w', markerfacecolor='darkblue', markersize=8, label='TPSA'),
    plt.Line2D([], [], marker='o', color='w', markerfacecolor='darkgreen', markersize=8, label='Num Rotatable Bonds'),
]
plt.legend(handles=legend_elements, loc='best', prop={'size': 10})

plt.show()
