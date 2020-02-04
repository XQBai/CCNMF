## Simulate the paired scRNA-seq and scDNA-seq datasets based on the common Clonal Genetic Copy
import numpy as np
import random
import pandas as pd

def Simulate_RNA(Configure):

    # The blow parameters are extracted from Splatter package "https://github.com/Irrationone/splatter" by
    # input the real paired  a triple-negative breast cancer xenograft scRNA-seq and scDNA-seq data from
    # <clonealign: statistical integration of independent single-cell RNA and DNA sequencing data from human cancers>.
    # where the data can be downloaded from "https://zenodo.org/record/2363826#.XjeDPKNVZE4"
    gene_mean_shape = 2.61
    gene_mean_rate = 0.97
    lib_loc = 8.6
    lib_scale = 0.592
    out_loc = 2.31
    out_scale = 0.743
    nGenes = Configure['geneticCN'][0].shape[0]
    nCells = Configure['geneticCN'][0].shape[1]
    config = {'RNA': []}
    for i in range(len(Configure)):
        outprob = Configure['Outlier'][i]
        dropout = Configure['Dropout'][i]
        # Simulated_RNA = np.zeros((nGenes, nCells))
        geneticCN = Configure['geneticCN'][i]
        # Gamma distribution
        gene_mean = np.random.gamma(gene_mean_shape, gene_mean_rate, size=(nGenes, 1))
        gene_mean = gene_mean/gene_mean.sum()
        One = np.ones((1, nCells))
        gene_mean = np.dot(gene_mean, One)
        # Dosage effect
        gene_mean_matrix = gene_mean * geneticCN
        # Library size
        L = np.random.lognormal(lib_loc, lib_scale, size=(1, nCells)) / 1152 * nCells
        #sum = np.sum(gene_mean[:, 0:10]) / 10
        #sum = np.sum(gene_mean_matrix[:, 0])
        L = np.ones((nGenes, 1)) * L
        Cell_gene_matrix = L * gene_mean_matrix

        Simulated_RNA = np.random.poisson(Cell_gene_matrix, size=(nGenes, nCells))


        Outlier_indictor = np.random.binomial(1, outprob, size=(nGenes, nCells))
        #Psi = np.random.lognormal(out_loc, out_scale, size = (1, nGenes))
        Simulated_RNA[np.where(Outlier_indictor == 1)] = Simulated_RNA[np.where(Outlier_indictor == 1)] + 2*np.round(np.mean(Simulated_RNA))
        # Dropout
        Dropout_indictor = np.random.binomial(1, dropout, size=(nGenes, nCells))
        Simulated_RNA[np.where(Dropout_indictor == 1)] = 0
        config['RNA'].append(Simulated_RNA)
    config = pd.Series(config['RNA'], name='RNA')
    return (config)


def Simulate_DNA(ProbMatrix, Configure):
    '''
    :param params
    :param ProbMatrix = [[0.42, 0.5, 0.08, 0, 0], [0.02, 0.52, 0.46, 0, 0],
    [0, 0, 0.5, 0.5, 0], [0, 0, 0.01, 0.4, 0.59], [0, 0, 0, 0.01, 0.99]]
    '''

    nGenes = Configure['geneticCN'][0].shape[0]
    nCells = Configure['geneticCN'][0].shape[1]

    index0 = np.array([0, 1, 2, 3, 4])
    index = index0.reshape(5, 1)
    config = {'CNV': []}
    for i in range(len(Configure)):
        outprob = Configure['Outlier'][i]
        dropout = Configure['Dropout'][i]
        Simulated_CNV = np.zeros((nGenes, nCells))
        geneticCN = Configure['geneticCN'][i]
        for m in range(nGenes):
            for n in range(nCells):
                if geneticCN[m, n] == 0:
                    a = np.random.multinomial(1, ProbMatrix[0], 1)
                    Simulated_CNV[m, n] = np.dot(a, index)[0][0]
                elif geneticCN[m, n] == 1:
                    a = np.random.multinomial(1, ProbMatrix[1], 1)
                    Simulated_CNV[m, n] = np.dot(a, index)[0][0]
                elif geneticCN[m, n] == 2:
                    a = np.random.multinomial(1, ProbMatrix[2], 1)
                    Simulated_CNV[m, n] = np.dot(a, index)[0][0]
                elif geneticCN[m, n] == 3:
                    a = np.random.multinomial(1, ProbMatrix[3], 1)
                    Simulated_CNV[m, n] = np.dot(a, index)[0][0]
                else:
                    a = np.random.multinomial(1, ProbMatrix[4], 1)
                    Simulated_CNV[m, n] = np.dot(a, index)[0][0]

        # utlier
        Outlier_indictor = np.random.binomial(1, outprob, size=(nGenes, nCells))
        Simulated_CNV[np.where(Outlier_indictor == 1)] = (Simulated_CNV +  np.round(2 * np.std(Simulated_CNV)))[np.where(Outlier_indictor == 1)]
        # Dropout
        Dropout_indictor = np.random.binomial(1, dropout, size=(nGenes, nCells))
        Simulated_CNV[np.where(Dropout_indictor == 1)] = 0
        config['CNV'].append(Simulated_CNV)
    config = pd.Series(config['CNV'], name='CNV')

    return (config)


def GeneticCN(d, nGenes, nCells):
    '''
    :param d: is a dict file. which is constructed by
    d = {'Ncluster' : [2, 3], 'Topology' : ['linear',  'bifurcate'],
        'C1Percent' : [0.5, 0.5] , 'C2Percent':[0.2, 0.4, 0.4],
        'Percentage' : [0.1, 0.2, 0.3]}
    :param nGenes: The number of genes
    :return:
    '''

    Ncluster = d['Ncluster']
    Topology = d['Topology']
    Percentage = d['Percentage']
    C1Percent = d['C1Percent']
    C2Percent = d['C2Percent']
    Outlier = d['Outlier']
    Dropout = d['Dropout']

    config = {'Ncluster': [], 'Topology': [], 'Percentage': [], 'Outlier':[], 'Dropout':[], 'V': [], 'geneticCN': []}
    for outlier in Outlier:
        for dropout in Dropout:
            for p in Percentage:
                for n in Ncluster:
                    V = np.zeros((np.max(Ncluster), nGenes))
                    if n == 2:
                        geneticCN = np.zeros((nGenes, nCells))
                        for i in range(n):

                            if i == 0:
                                V[i, :] = 2
                                x = int(nCells * C1Percent[0])
                                V1 = V[i, :].reshape(len(V[i, :]), 1)
                                One = np.ones(x)
                                One = One.reshape(1, x)
                                geneticCN[:, 0:x] = np.dot(V1, One)
                            elif i == 1:
                                loc = random.sample(list(range(nGenes)), int(nGenes * p))
                                V[i, :] = 2
                                V[i, loc] = np.random.choice([0, 1, 3, 4], len(loc), replace=True)
                                x1 = int(nCells * C1Percent[0])
                                x2 = int(nCells * C1Percent[1])
                                V1 = V[i, :].reshape(len(V[i, :]), 1)
                                One = np.ones(x2)
                                One = One.reshape(1, x2)
                                geneticCN[:, x1: (x1 + x2)] = np.dot(V1, One)
                        config['V'].append(V[0:i + 1, :])
                        config['Percentage'].append(p)
                        config['Ncluster'].append(n)
                        config['Topology'].append('NA')
                        config['geneticCN'].append(geneticCN)

                    elif n == 3:
                        for t in Topology:
                            V = np.zeros((np.max(Ncluster), nGenes))
                            geneticCN = np.zeros((nGenes, nCells))
                            for i in range(n):
                                # geneticCN = np.zeros((nGenes, nCells))
                                if i == 0:
                                    V[i, :] = 2
                                    x = int(nCells * C2Percent[0])
                                    V1 = V[i, :].reshape(len(V[i, :]), 1)
                                    One = np.ones(x)
                                    One = One.reshape(1, x)
                                    geneticCN[:, 0:x] = np.dot(V1, One)
                                elif i == 1:
                                    loc = random.sample(list(range(nGenes)), int(nGenes * p))
                                    V[i, :] = 2
                                    V[i, loc] = np.random.choice([0, 1, 3, 4], len(loc), replace=True)
                                    x1 = int(nCells * C2Percent[0])
                                    x2 = int(nCells * C2Percent[1])
                                    V1 = V[i, :].reshape(len(V[i, :]), 1)
                                    One = np.ones(x2)
                                    One = One.reshape(1, x2)
                                    geneticCN[:, x1: (x1 + x2)] = np.dot(V1, One)
                                else:
                                    if t == 'linear':
                                        V[i, :] = V[i - 1, :]
                                        loc = random.sample(list(np.where(V[i, :] == 2)[0]), int(nGenes * p))
                                        V[i, loc] = np.random.choice([0, 1, 3, 4], len(loc), replace=True)
                                        x3 = int(nCells * C2Percent[2])
                                        V1 = V[i, :].reshape(len(V[i, :]), 1)
                                        One = np.ones(x3)
                                        One = One.reshape(1, x3)
                                        geneticCN[:, (nCells - x3):] = np.dot(V1, One)
                                    elif t == 'bifurcate':
                                        V[i, :] = V[i - 2, :]
                                        loc = random.sample(list(np.where(V[i - 1, :] == 2)[0]), int(nGenes * p))
                                        V[i, loc] = np.random.choice([0, 1, 3, 4], len(loc), replace=True)
                                        x3 = int(nCells * C2Percent[2])
                                        V1 = V[i, :].reshape(len(V[i, :]), 1)
                                        One = np.ones(x3)
                                        One = One.reshape(1, x3)
                                        geneticCN[:, (nCells - x3):] = np.dot(V1, One)
                            config['V'].append(V[0:i + 1, :])
                            config['Percentage'].append(p)
                            config['Ncluster'].append(n)
                            config['Topology'].append(t)
                            config['Outlier'].append(outlier)
                            config['Dropout'].append(dropout)
                            config['geneticCN'].append(geneticCN)

    Configure = pd.DataFrame(config, index=['ID' + str(i + 1) for i in range(len(config['Ncluster']))])
    return (Configure)
