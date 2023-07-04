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
for i, smi in enumerate(smiles):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        print(f'Invalid SMILES string: {smi}')
    else:
        mols.append(mol)
        valid_logP.append(logP[i])

# Calculate molecular descriptors
X = np.array([calc_descriptors(mol) for mol in mols])

# Perform PCA to reduce the dimensionality of the data
pca = PCA(n_components=3)
X_pca = pca.fit_transform(X)

# Create a 3D plot with a black background and colored points
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
sc = ax.scatter(X_pca[:, 0], X_pca[:, 1], X_pca[:, 2], c=valid_logP, cmap='viridis')
ax.set_facecolor('black')
ax.set_xlabel('PC1')
ax.set_ylabel('PC2')
ax.set_zlabel('PC3')
plt.colorbar(sc, label='LogP')
plt.title('Chemical Space')

# Add a legend to show which molecular descriptors were used
legend_elements = [
    plt.Line2D([], [], marker='', linestyle='', label='Molecular Weight'),
    plt.Line2D([], [], marker='', linestyle='', label='TPSA'),
    plt.Line2D([], [], marker='', linestyle='', label='Num Rotatable Bonds'),
]
plt.legend(handles=legend_elements, loc='best')

plt.show()
