import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense
import argparse

parser = argparse.ArgumentParser(description='ML Model')
parser.add_argument('-i', type=str, default=None, required=True, help='File input (JSON)')
args = parser.parse_args()

primes = {
    'A': 2,
    'C': 3,
    'G': 5,
    'T': 7,
    'R': 2 * 5,
    'Y': 3 * 7,
    'S': 5 * 3,
    'W': 2 * 7,
    'K': 5 * 7,
    'M': 2 * 3,
    'B': 3 * 5 * 7,
    'D': 2 * 5 * 7,
    'H': 2 * 3 * 7,
    'V': 2 * 3 * 5,
    'N': 2 * 3 * 5 * 7
}

weights = {
    'A': 4,
    'C': 4,
    'G': 4,
    'T': 4,
    'R': 3,
    'Y': 3,
    'S': 3,
    'W': 3,
    'K': 3,
    'M': 3,
    'B': 2,
    'D': 2,
    'H': 2,
    'V': 2,
    'N': 1,
    '-': 1
}

# Assume you have a DataFrame `df` with your data
df = pd.read_csv(args.i)

# Prepare the dataset
X = df[['m_score', 'r_score', 'smith_waterman', 'similarity', 'enzymes_nearby']].values
y = df['target'].values

# Normalize features
scaler = StandardScaler()
X = scaler.fit_transform(X)

# Split into train and test sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Define the model
model = Sequential()
model.add(Dense(10, input_dim=X_train.shape[1], activation='relu'))  # Input layer
model.add(Dense(10, activation='relu'))  # Hidden layer
model.add(Dense(1, activation='sigmoid'))  # Output layer

# Compile the model
model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['accuracy'])

# Train the model
model.fit(X_train, y_train, epochs=10, batch_size=32, verbose=1)

# Evaluate the model
loss, accuracy = model.evaluate(X_test, y_test, verbose=1)
print(f'Test Accuracy: {accuracy * 100:.2f}%')
