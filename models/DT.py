import pandas as pd
from sklearn import tree
import sys


# argv1 - training set in feather format
# argv2 - test set in feather format
# argv3 - output in feather format

X, y = ["coverage", "fetal_fraction", "euploidy", "trisomy"], "condition"  # 2-state model
# X, y = ["coverage", "fetal_fraction", "euploidy", "trisomy", "paternal_trisomy"], "genotype"  # 7-state model

# fit model
training = pd.read_feather(sys.argv[1])
training_X = training[X]
training_y = training[y]
clf = tree.DecisionTreeClassifier(max_depth = 3, random_state = 123)
clf = clf.fit(training_X, training_y)

# predict
test = pd.read_feather(sys.argv[2])
test_X = test[X]
test_y = clf.predict(test_X)

# write results to file
results = test.assign(prediction = test_y)
results.to_feather(sys.argv[3])
