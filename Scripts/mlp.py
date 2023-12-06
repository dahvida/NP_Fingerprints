"""
Implementation of a Multilayer Perceptron in Pytorch so that it is compatible
with sklearn's API
"""

import torch
from torch import nn
import numpy as np
from torch.utils.data import Dataset, DataLoader

class MLP:

    def __init__(self,
                 input_size,
                 hidden_size,
                 dropout_rate,
                 learning_rate,
                 batch_size 
                 ):

        # store input in dict
        self.param_dict = {
            "input_size": input_size,                               # fingerprint dimensionality 
            "hidden_size": hidden_size,                             # number of neurons per layer
            "dropout": dropout_rate,                                # dropout rate for LSTM encoders
            "batch_size": batch_size,                               # samples per gradient step
            "learning_rate": learning_rate,                         # update step magnitude
            "epochs": 100
        }

        # make model architecture
        MLP_list = [] 
        MLP_list.append(nn.Linear(input_size, hidden_size))
        MLP_list.append(nn.BatchNorm1d(hidden_size))
        MLP_list.append(nn.Dropout(dropout_rate))
        MLP_list.append(nn.ReLU())
        MLP_list.append(nn.Linear(hidden_size, hidden_size))
        MLP_list.append(nn.BatchNorm1d(hidden_size))
        MLP_list.append(nn.Dropout(dropout_rate))
        MLP_list.append(nn.ReLU())
        MLP_list.append(nn.Linear(hidden_size, 1))
        MLP_list.append(nn.Sigmoid())
        self.MLP = nn.Sequential(*MLP_list)

        # make optimizer
        self.optimizer = torch.optim.AdamW(
            list(self.MLP.parameters()),
            lr=learning_rate,
            )        
    
    def forward(self, x):
        return self.MLP(x)

    def fit(self, x, y):
        
        self.threshold = np.mean(y)

        self.MLP.train()
        loss = nn.BCELoss()
        train_data = DataGenerator(x, y)
        train_data = DataLoader(train_data, 
                                batch_size = self.param_dict["batch_size"],
                                shuffle = True)
        epochs = self.param_dict["epochs"]
        
        for epoch in range(epochs):
            
            # initialize loss containers
            train_loss = 0.0
            
            # loop over training set
            for batch_idx, (x_batch, y_batch) in enumerate(train_data):
                
                # reset grad
                self.optimizer.zero_grad()

                # get predictions
                preds = self.forward(x_batch)[:,0]
                
                # get loss
                batch = loss(preds, y_batch)
                
                # do backwards step
                batch.backward()
                self.optimizer.step()
                
                # add i-th loss to training container
                train_loss += batch.item()

            # Print mean epoch loss
            print(f"Epoch [{epoch + 1}/{epochs}], Train loss: {train_loss / len(train_data):.3f}")

    def predict(self, x):

        self.MLP.eval()
        x = torch.tensor(x, dtype=torch.float32)
        preds = self.forward(x)
        preds = preds.detach().numpy()
        preds = (preds >= self.threshold).astype(int)

        return preds

    def predict_proba(self, x):
        self.MLP.eval()
        x = torch.tensor(x, dtype=torch.float32)
        preds = self.forward(x)
        preds = preds.detach().numpy()
        preds = np.concatenate((1-preds, preds), axis=1)

        return preds 
        
class DataGenerator(Dataset):

    def __init__(self,
                 x,
                 y,
                 ):
        self.x = x
        self.y = y

    def __len__(self):
        """
        Method necessary for Pytorch training
        """
        return len(self.x)

    def __getitem__(self, idx):
        """
        Method necessary for Pytorch training
        """
        
        x_i = torch.tensor(self.x[idx],
                             dtype=torch.float32)
        y_i = torch.tensor(self.y[idx],
                             dtype=torch.float32)

        return x_i, y_i


