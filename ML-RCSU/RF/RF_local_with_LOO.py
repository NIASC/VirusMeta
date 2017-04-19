# Databricks notebook source exported at Wed, 21 Dec 2016 12:41:17 UTC

# you should have sklearn 0.18!!

import sklearn

print sklearn.__version__

# http://matplotlib.org/faq/howto_faq.html#generate-images-without-having-a-window-appear
import matplotlib as mpl

mpl.use('Agg')
import matplotlib.pyplot as plt

# COMMAND ----------

# MAGIC %md
# MAGIC # Split into training and test set
# MAGIC
# MAGIC Different helper functions for splitting

# COMMAND ----------

import numpy as np
from numpy import loadtxt
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score

# the normal random shuffle split
def TT_split(features, labels, prop=0.9):
    nr_samples = features.shape[0]
    nr_train = int(nr_samples * prop)
    indx = range(nr_samples)
    np.random.shuffle(indx)
    features_temp = features[indx, :]
    labels_temp = labels[indx]

    train_X = features_temp[:nr_train, :]
    train_Y = labels_temp[:nr_train]
    test_X = features_temp[nr_train:, :]
    test_Y = labels_temp[nr_train:]
    return train_X, train_Y, test_X, test_Y


# split by giving an array of experiment numbers and the nr to be left to the test set
def project_split(features, labels, projects, experiment_nr):
    projects = np.array(projects)
    print "test set size", np.shape(np.where(projects == experiment_nr)[0])
    assert (len(np.where(projects != experiment_nr)[0]) > 0)

    train_X = features[np.where(projects != experiment_nr)[0]]
    train_Y = labels[np.where(projects != experiment_nr)[0]]
    test_X = features[np.where(projects == experiment_nr)[0]]
    test_Y = labels[np.where(projects == experiment_nr)[0]]
    return train_X, train_Y, test_X, test_Y


# helper to get the project ID from the sequenceID (which is a string)
def seqID_to_projects(seqID):
    projects = np.array([x[x.find("rr"):] for x in seqID])
    # print np.unique(projects), np.shape(np.unique(projects))
    exp_names = np.unique(projects)
    pr_nr = np.zeros(np.shape(projects))
    for i, name in enumerate(exp_names):
        print i, name, np.shape(np.where(projects == name))
        pr_nr[np.where(projects == name)[0]] = i
    return pr_nr


# COMMAND ----------

# MAGIC %md # # Simple 90/10 split

# COMMAND ----------

from numpy import loadtxt
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score

f = ["data/after_clust_concat_al2_"]  # list of datasets to use

forests = []
for i, folder in enumerate(f):
    features = loadtxt(folder + "matrix.txt")
    labels = loadtxt(folder + "values.txt", dtype=int)
    print features.shape, labels.shape

    [train_X, train_Y, test_X, test_Y] = TT_split(features, labels)
    print "shapes after split ", train_X.shape, train_Y.shape, test_X.shape, test_Y.shape
    print "class counts in test set: ", np.bincount(test_Y)

    # rf = RandomForestClassifier(class_weight="balanced")
    rf = RandomForestClassifier(n_estimators=100)
    rf.fit(train_X, train_Y)

    # print "feature imp:",rf.feature_importances_
    pred = rf.predict(test_X)
    prob = rf.predict_proba(test_X)
    prs = []
    res = []
    for th in np.arange(0.5, 0.991, 0.05):
        # print "\r",th
        not_sure = np.where((np.max(prob, axis=1) < th))[0]
        pred_th = pred.copy()
        pred_th[not_sure] = 2
        pr = precision_score(np.concatenate((test_Y, [2])), np.concatenate((pred_th, [2])), average=None)[1]
        re = recall_score(np.concatenate((test_Y, [2])), np.concatenate((pred_th, [2])), average=None)[1]
        print precision_score(np.concatenate((test_Y, [2])), np.concatenate((pred_th, [2])), average=None)
        print recall_score(np.concatenate((test_Y, [2])), np.concatenate((pred_th, [2])), average=None)
        prs.append(pr)
        res.append(re)
    report = classification_report(test_Y, pred, target_names=["not virus", "virus", "not_classified"])
    print "report at 0.5: ", report
    # report99 = classification_report(test_Y, pred_th, target_names=["not virus", "virus","not_classified"])
    # print "report at 0.99: ",  report99
    plt.plot(np.array(res) * 100, 100 * np.array(prs))
    print "----------------------------------------------"
    print np.shape(prs)
    plt.ylabel("precision")
    plt.axhline(90, color="r", linestyle="dashed")
    plt.xlabel("recall")
    # plt.ylim([85, 101])
    # plt.xlim([10, 45])
    # plt.savefig("plot_nr"+str(i)+".png")
    plt.show()
    plt.clf()
    forests.append(rf)

# COMMAND ----------

# LEAVE ONE OUT

f = ["data/after_clust_concat_al2_"]

# save the forests and their statistics in these variables
precisions_for_LOO = [[], []]
recalls_for_LOO = [[], []]
class_counts = [[], []]

thresholds = np.arange(0.51, 1.0, 0.01)  # THIS SAME PARAM will be used in the next cell for plotting as well
len_thresholds = len(thresholds)

# for different datasets
for i, folder in enumerate(f):
    print "---------#################  LEAVE ONE OUT  ###################--------\n", folder
    features = loadtxt(folder + "matrix.txt")
    labels = loadtxt(folder + "values.txt", dtype=int)
    projects = loadtxt(folder + "projects.txt", dtype=int)

    print "overall data shape: ", features.shape, labels.shape, projects.shape, "\n overall counts of classes ", np.bincount(
        labels)

    nr_of_projects = len(np.unique(projects))
    print "We have ", nr_of_projects, " different projects in the dataset"

    # for each project we train a new model
    for LO in range(nr_of_projects):
        [train_X, train_Y, test_X, test_Y] = project_split(features, labels, projects, LO)

        print "LOO " + str(LO) + " shapes after split ", train_X.shape, train_Y.shape, test_X.shape, test_Y.shape
        print "counts in test set", np.bincount(test_Y)
        cc = np.bincount(test_Y)

        if len(cc) == 1:  # if there is only one label to be counted
            print cc, "there are NO VIRUSES in the test set!! (no training, instead adding 0-s by hand)"
            class_counts[i].append(np.array([cc[0], 0]))
            precisions_for_LOO[i].append([0] * len_thresholds)
            recalls_for_LOO[i].append([0] * len_thresholds)
            continue  # end this iteration in the FOR-loop

        class_counts[i].append(cc)
        rf = RandomForestClassifier(n_estimators=100)
        rf.fit(train_X, train_Y)

        ## from sklearn.externals import joblib
        ## joblib.dump(clf, 'filename.pkl')
        ## #then your colleagues can load it
        ## clf = joblib.load('filename.pk1')

        ## #OR:
        ## import cPickle
        ## rf = RandomForestRegresor()
        ## rf.fit(X, y)
        ## with open('path/to/file', 'wb') as f:
        ##     cPickle.dump(rf, f)
        ## # in your prediction file
        ## with open('path/to/file', 'rb') as f:
        ##     rf = cPickle.load(f)
        ## preds = rf.predict(new_X)

        ## #OR:
        ## import dill
        ## wd = "/whatever/you/want/your/working/directory/to/be/"
        ## rf= RandomForestRegressor(n_estimators=250, max_features=9,compute_importances=True)
        ## rf.fit(Predx, Predy)
        ## dill.dump(rf, open(wd + "filename.obj","wb"))
        ## with open(wd + "filename.obj","wb") as f:
        ##     dill.dump(rf,f)
        ## model = dill.load(open(wd + "filename.obj","rb"))

        pred = rf.predict(test_X)
        prob = rf.predict_proba(test_X)
        prs = []
        res = []
        for th in thresholds:
            # print "\r",th
            not_sure = np.where((np.max(prob, axis=1) < th))[0]
            pred_th = pred.copy()
            pred_th[not_sure] = 2
            pr = precision_score(test_Y, pred_th, average=None)
            re = recall_score(test_Y, pred_th, average=None)
            if len(pr) == 1:
                prs.append(0)
                res.append(0)
            else:
                prs.append(pr[1])
                res.append(re[1])
        precisions_for_LOO[i].append(prs)
        recalls_for_LOO[i].append(res)

        # report = classification_report(test_Y, pred, target_names=["not virus", "virus","not_classified"])
        # print "report at 0.5: ", report
        # report99 = classification_report(test_Y, pred_th, target_names=["not virus", "virus","not_classified"])
        # print "report at 0.99: ",  report99

# COMMAND ----------

# visualize the statistics (read from a saved pickle)


# TP+FN=count_vir
# RE=TP/(TP+FN) ---> TP = RE*count_vir && FN= count_vir-TP
# PR=TP/(TP+FP) ---> FP = TP/PR -TP

# thresholds = np.arange(0.51,1.0,0.01)

plt.clf()
prcs = []
recs = []
counts = np.array(class_counts[0])

print counts, counts.shape, counts[1, 1]

for th in range(len(thresholds)):  # seven thresholds
    sum_TP = 0
    sum_FP = 0
    sum_FN = 0
    sum_count_vir = 0
    for loo in range(len(counts)):  # 20 experimetns?
        count_vir = counts[loo, 1]
        # print count_vir
        sum_count_vir += count_vir

        RE = recalls_for_LOO[0][loo][th]
        PR = precisions_for_LOO[0][loo][th]
        TP = RE * count_vir
        FN = count_vir - TP
        if TP == 0 or PR == 0:  # experiments with so few samples that there are no predicted viruses
            # print "have to skip step ", loo, "where count is ", count_vir
            continue  # do not consider this (do not add to sum)
        else:
            FP = TP / PR - TP
            sum_TP += TP
            sum_FP += FP
            sum_FN += FN
            # print "count", count_vir,"RE", RE,"PR",PR,"TP", TP,"FP",FP,"FN",FN

    if not sum_FN == 0:
        print "for threshold ", thresholds[
            th], " sum_TP, sum_FP, sum_FN, sum_count_vir ", sum_TP, sum_FP, sum_FN, sum_count_vir
        print "resulting in PR, RE: ", sum_TP / (sum_TP + sum_FP), sum_TP / sum_count_vir
        prcs.append(sum_TP / (sum_TP + sum_FP))
        recs.append(sum_TP / sum_count_vir)
        print "#############################"
    else:
        prcs.append(0)
        recs.append(0)

# l,= plt.plot(prcs,recs)
plt.plot(100 * np.array(recs), 100 * np.array(prcs))
plt.scatter(100 * np.array(recs), 100 * np.array(prcs))
plt.ylim((62, 100))
plt.ylabel("Precision")
plt.xlabel("Recall")
plt.title("Precision-Recall trade-off curve in leave-out-experiments")
plt.axvline(26, color="r", linestyle="dashed")
plt.axhline(90, color="r", linestyle="dashed")
plt.show()
plt.clf()

# COMMAND ----------

# In this plot the numbers are written by hand as of now...

plt.clf()
# write results here:
plt.bar([2, 1], [93, 85], 0.3)  # precision of nonvirus and virus
plt.bar([2.35, 1.35], [99, 47], 0.3, color="green")  # recall of nonvirus and virus

# no need to change the following
plt.ylabel("Precision/Recall in %")
plt.xticks([2.325, 1.325], ["Non-virus", "Virus"])
plt.xlim([0.95, 2.75])
plt.ylim([0, 118])
plt.legend(["Precision", "Recall"], loc=9)
plt.axhline(100, color="gray", linestyle="dashed")
plt.title("Precision and recall for 90/10 split with RF")
plt.show()

# COMMAND ----------

# RESULTS FROM NNETWORK, Hardcoded, do not touch!

plt.clf()
plt.bar([2, 1], [94, 78], 0.3)  # precision
plt.bar([2.35, 1.35], [98, 52], 0.3, color="green")  # recall

plt.ylabel("Precision/Recall in %")
plt.xticks([2.325, 1.325], ["Non-virus", "Virus"])
plt.xlim([0.95, 2.75])
plt.ylim([0, 118])
plt.legend(["Precision", "Recall"], loc=9)
plt.axhline(100, color="gray", linestyle="dashed")
plt.title("Precision and recall for 90/10 split with feed-forward NN")
plt.show()

# COMMAND ----------

# MAGIC %md # TRAIN ON DNA, TEST ON SOMETHING ELSE
# MAGIC
# MAGIC
# MAGIC NOTICE: In here you can replace the "TEST" with any other dataset

# COMMAND ----------

DNA = "data/after_clust_concat_al2_"
TEST = "data/after_clust_concat_al2_"  ## LEFT IT THE SAME AS TEST, replace it with desired dataset

DNA_features = loadtxt(folder + "matrix.txt")
DNA_labels = loadtxt(folder + "values.txt", dtype=int)
print DNA_features.shape, DNA_labels.shape
###############

TEST_features = loadtxt(folder + "matrix.txt")
TEST_labels = loadtxt(folder + "values.txt", dtype=int)
print TEST_features.shape, TEST_labels.shape

# Train model on whole DNA dataset, then test on TEST
rf = RandomForestClassifier(n_estimators=100)
rf.fit(DNA_features, DNA_labels)

importance = rf.feature_importances_
print importance

# test the models on some dataset
pred = rf.predict(TEST_features)
prob = rf.predict_proba(TEST_features)

report = classification_report(TEST_labels, pred, target_names=["not virus", "virus"])
print "Report on test set at threshold 0.5: \n", report
print "#################################################"

# COMMAND ----------

# MAGIC %md ## JUST LITTLE THINGS TO PLOT

# COMMAND ----------

# not virus       0.93      0.98      0.95     19132
#     virus       0.75      0.49      0.59      2791
plt.clf()

# RF
plt.bar([1, 1.8], [93, 85], 0.1)
plt.bar([1.11, 1.91], [99, 48], 0.1, color="green")

# NNET
plt.bar([1.25, 2.05], [94, 78], 0.1, color="lightskyblue")
plt.bar([1.36, 2.16], [98, 52], 0.1, color="darksage")

plt.ylabel("Precision/Recall in %")
plt.xticks([1.225, 2.025], ["Non-virus", "Virus"])
plt.xlim([0.95, 2.8])
plt.ylim([0, 118])
plt.legend(["Precision RF", "Recall RF", "Precision NN", "Recall NN"])
plt.axhline(100, color="gray", linestyle="dashed")
plt.show()

# COMMAND ----------

# these are feature importances for clustered_concatenated_2ORF
feature_imp = [
    [0.00987578, 0.01074493, 0.01267057, 0.01142763, 0.01653067, 0.03889236,
     0.01185394, 0.01216873, 0.01229994, 0.01298374, 0.01398363, 0.01297234,
     0.01246973, 0.01620767, 0.02466992, 0.01760757, 0.0127776, 0.01411569,
     0.013858, 0.03993107, 0.01175226, 0.01169014, 0.02314367, 0.01574593,
     0.02629938, 0.01408819, 0.01263433, 0.01306898, 0.01341817, 0.01617842,
     0.02420086, 0.01332598, 0.01083113, 0.01095561, 0.01899907, 0.01673964,
     0.01673531, 0.01522593, 0.02040868, 0.02329504, 0.0156791, 0.01531973,
     0.01280192, 0.0144967, 0.01492245, 0.0320429, 0.03227925, 0.01534282,
     0.02538328, 0.02622538, 0.01953664, 0.02318912, 0.01098458, 0.01414962,
     0.01436501, 0.01215648, 0.01676019, 0.01187635, 0.01571022]]

## this defines a color scheme
tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]

# Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.
for i in range(len(tableau20)):
    r, g, b = tableau20[i]
    tableau20[i] = (r / 255., g / 255., b / 255.)

tableau20 = np.array(tableau20)

codons = ["Cys", "Cys", "Asp", "Asp", "Ser", "Ser", "Ser", "Ser", "Ser", "Ser", "Gln", "Gln", "Lys", "Lys",
          "Gly", "Gly", "Gly", "Gly", "Pro", "Pro", "Pro", "Pro", "Thr", "Thr", "Thr", "Thr", "Phe", "Phe",
          "Ala", "Ala", "Ala", "Ala", "His", "His", "Ile", "Ile", "Ile", "Glu", "Glu",
          "Leu", "Leu", "Leu", "Leu", "Leu", "Leu", "Arg", "Arg", "Arg", "Arg", "Arg", "Arg",
          "Val", "Val", "Val", "Val", "Asn", "Asn", "Tyr", "Tyr"]

temporary = "Cys"
aid = 0
amino_id = []
aminos = ["Cys"]
for i in range(59):
    if codons[i] == temporary:
        amino_id.append(aid)
    else:
        aid += 1
        temporary = codons[i]
        amino_id.append(aid)
        aminos.append(codons[i])
print amino_id

plt.scatter((range(59)), feature_imp[-1], color=tableau20[amino_id, :])
plt.legend(aminos)
plt.xlim([-1, 60])

import matplotlib.patches as mpatches

recs = []
for i in range(0, len(aminos)):
    recs.append(mpatches.Rectangle((0, 0), 1, 1, fc=tableau20[i, :]))
plt.legend(recs, aminos, loc=9, ncol=6, fontsize=8)
plt.title("Importance of features in the Random Forest classifier")
plt.show()

# COMMAND ----------


