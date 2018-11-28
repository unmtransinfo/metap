
# Use of XGBoost
(XGBoost = eXtreme Gradient Boosting)
<img src="https://raw.githubusercontent.com/dmlc/dmlc.github.io/master/img/logo-m/xgboost.png" width="135" style="float:right" />
 * <https://xgboost.ai/>
 * [GitHub](https://github.com/dmlc/xgboost)
 * [Documentation](https://xgboost.readthedocs.io/en/latest/)
 * [R package](https://cran.r-project.org/web/packages/xgboost/)

MetapathML employs XGBoost via the R package API. The inputs to XGBoost are datasets specific to each disease or phenotype. For each disease/phenotype some known associated genes correspond with the positive Y labels in the dataset. XGBoost parameters are optimized via grid search, i.e. iterative testing over discrete parameter value combinations. XGBoost is used in tree booster mode, where the following parameters are tuned:

parameter | choices | default | description
---: | :---: | :---: | :---
max_depth | (5, 7, 10) | 6 | Maximum depth of a tree. Increasing this value will make the model more complex and more likely to overfit.
eta | (.05, .1, .15, .20) | .3 | Step size shrinkage used in update to prevent overfitting.
gamma | (.01, .1, 1) | 0 | Minimum loss reduction required to make a further partition on a leaf node of the tree. The larger gamma is, the more conservative the algorithm will be.
min_child_weight | (0, 1, 2) | 1 | Minimum sum of instance weight (hessian) needed in a child.
subsample | (.8, .9, 1) | 1 | Subsample ratio of the training instances. Setting it to 0.5 means that XGBoost would randomly sample half of the training data prior to growing trees. and this will prevent overfitting.
colsample_bytree | (.5, .6, .7, .8, .9, 1) | 1 | Subsample ratio of columns when constructing each tree.

Note that several parameters are not optimized. Some listed here, with focus on overfitting/regularization.

parameter | value | default | description
---: | :---: | :---: | :---
objective | binary:logistic | reg:linear | Specify the learning task and the corresponding learning objective.
metrics | auc |  | Evaluation metric[s] for validation data.
max_delta_step | | 0 | Maximum delta step we allow each leaf output to be. If the value is set to 0, it means there is no constraint.
colsample_bylevel | | 1 | Subsample ratio of columns for each split, in each level.
lambda | | 1 | L2 regularization term on weights. Increasing this value will make model more conservative.
alpha | | 0 | L1 regularization term on weights. Increasing this value will make model more conservative.
tree_method | | auto | The tree construction algorithm used in XGBoost. See description in the [reference paper](http://arxiv.org/abs/1603.02754).
scale_pos_weight | sumneg/sumpos | 1 | Control the balance of positive and negative weights, useful for unbalanced classes. A typical value to consider: sum(negative instances) / sum(positive instances). See [Parameters Tuning](https://xgboost.readthedocs.io/en/latest/tutorials/param_tuning.html) for more discussion.


See <https://xgboost.readthedocs.io/en/latest/parameter.html> for parameter details.
