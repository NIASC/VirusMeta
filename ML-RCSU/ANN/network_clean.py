#!/usr/bin/env python

# https://www.csc.fi/web/training/-/yandex_2017
################ VERSIONS of packages I am using
# Python 2.7.10
# keras.__version__ '1.1.0'
# numpy.__version__ '1.11.1'
# sklearn.__version__ '0.18'
# matplotlib.__version__ '1.5.3'
# argparse.__version__ '1.1'

import numpy as np
from keras.datasets import mnist
from keras.models import Sequential
from keras.layers.core import Dense, Dropout, Activation
from keras.optimizers import SGD, Adam, RMSprop
from keras.callbacks import ModelCheckpoint, LearningRateScheduler, Callback
from keras.utils import np_utils
import keras.backend as K
from itertools import product
from functools import partial
from sklearn.cross_validation import train_test_split
import argparse
import sklearn.metrics as metrics
from sklearn.metrics import classification_report

# http://matplotlib.org/faq/howto_faq.html#generate-images-without-having-a-window-appear
import matplotlib as mpl

mpl.use('Agg')
import matplotlib.pyplot as plt

# this is a helper function that allows us to define any class weights we want (in a n x n matrix)
def w_categorical_crossentropy(y_true, y_pred, weights):
    nb_cl = len(weights)
    final_mask = K.zeros_like(y_pred[:, 0])
    y_pred_max = K.max(y_pred, axis=1)
    y_pred_max = K.expand_dims(y_pred_max, 1)
    y_pred_max_mat = K.equal(y_pred, y_pred_max)
    for c_p, c_t in product(range(nb_cl), range(nb_cl)):
        final_mask += (
        K.cast(weights[c_t, c_p], K.floatx()) * K.cast(y_pred_max_mat[:, c_p], K.floatx()) * K.cast(y_true[:, c_t],
                                                                                                    K.floatx()))
    return K.categorical_crossentropy(y_pred, y_true) * final_mask

# this is the main object - MetaGenomics FeedForward classifier
class MetaFF:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

    # CREATING THE NETWORK
    def init(self, nb_inputs, nb_outputs):
        print "Creating model..."

        self.model = Sequential()
        # create n hidden layers
        for i in xrange(self.layers):
            if i == 0:
                layer = Dense(self.hidden_nodes, activation=self.activation, input_shape=(nb_inputs,))
            else:
                layer = Dense(self.hidden_nodes, activation=self.activation)
            self.model.add(layer)
            if self.dropout > 0:  # dropout after each layer
                self.model.add(Dropout(self.dropout))
        # and the final layer
        self.model.add(Dense(nb_outputs))  # no nonlinearity here before softmax
        self.model.add(Activation('softmax'))
        self.model.summary()  # prints a summary

        if args.class_weight_power is not None:  # playing with class weights, using the custom loss function created above
            ncce = partial(w_categorical_crossentropy, weights=np.array(
                [[args.class_weight_power[0], args.class_weight_power[0]],
                 [args.class_weight_power[1], args.class_weight_power[1]]]))
            ncce.__name__ = 'w_categorical_crossentropy'
            self.model.compile(loss=ncce, optimizer=self.optimizer)
        else:  # if not using class weights
            print "Compiling model..."
            self.model.compile(loss='categorical_crossentropy', optimizer=self.optimizer)
            # end of defining the model

    # TRAINING
    def fit(self, train_X, train_y, valid_X, valid_y, save_path):
        # saving the model, normally (save_best_only=False) we do it after every epoch
        callbacks = [ModelCheckpoint(filepath=save_path, verbose=1, save_best_only=args.save_best_model_only)]

        if self.lr_epochs > 0:  # create LR scheduler if needed
            def lr_scheduler(epoch):
                lr = self.lr * self.lr_factor ** int(epoch / self.lr_epochs)
                print "Epoch %d: learning rate %g" % (epoch + 1, lr)
                return lr

            callbacks.append(LearningRateScheduler(lr_scheduler))

        history = self.model.fit(train_X, train_y, batch_size=self.batch_size, nb_epoch=self.epochs,
                                 validation_data=(valid_X, valid_y),
                                 shuffle=True, verbose=self.verbose, callbacks=callbacks)

    # EVALUATING the model
    def eval(self, X, y, load_path):
        self.model.load_weights(load_path)  # this works even when evaluating during training

        pred_y = self.model.predict(X, batch_size=self.batch_size)
        pred_y_labels = np.argmax(pred_y, axis=1)

        if not args.binary:
            report = classification_report(y, pred_y_labels,
                                           target_names=["bacteria", "inverte", "vertebrate", "plant", "virus"])
        else:
            report = classification_report(y, pred_y_labels, target_names=["not virus", "virus"])
        return report

    def predict(self, X, load_path):
        self.model.load_weights(load_path)
        pred_y = self.model.predict(X, batch_size=self.batch_size)
        print 'pred_y shape:', pred_y.shape
        return pred_y

    def predict_and_save(self, unknown, load_path, pp_result):
        self.model.load_weights(load_path)
        pred_y = self.model.predict(unknown, batch_size=self.batch_size)  # returns [P_non,P_vir]
        print 'plot threshold: pred_y shape ', pred_y.shape

        # we will create list with prediction and probability
        pred_and_prob = np.vstack([np.argmax(pred_y, axis=1), np.max(pred_y, axis=1)])
        # print np.min(np.max(pred_y,axis=1))
        pred_and_prob = np.transpose(pred_and_prob)
        # save pp
        with open(pp_result, "w") as self_pp_result:
            # self_pp_result.write(pred_and_prob)
            np.savetxt(self_pp_result, pred_and_prob)

    # THRESHOLDING by output (does not reflect uncertainty)
    def plot_threshold(self, X, y, load_path):
        self.model.load_weights(load_path)
        pred_y = self.model.predict(X, batch_size=self.batch_size)  # returns [P_non,P_vir]
        print 'plot threshold: pred_y shape ', pred_y.shape

        # we will create list with prediction and probability
        pred_and_prob = np.vstack([np.argmax(pred_y, axis=1), np.max(pred_y, axis=1)])
        # print np.min(np.max(pred_y,axis=1))
        pred_and_prob = np.transpose(pred_and_prob)

        pp = []
        for threshold in np.arange(0.5, 0.991, 0.01):  # for different thresholds
            pr_and_pr = pred_and_prob.copy()
            for i in range(pred_and_prob.shape[0]):  # for all datapoints
                if pr_and_pr[i, 1] < threshold:  # if proba below threshold
                    pr_and_pr[i, 0] = -1  # set prediction to "no prediction"

            above = np.where(pr_and_pr[:, 0] > -1)[0]  # find the indexes of all instances above thresh
            trues = y[above]
            preds = pr_and_pr[above, 0].flatten()

            pr = metrics.precision_score(trues, preds, average=None)  # returns precision for each class (2 values)
            tr = (trues == preds)
            tp = trues * tr

            re = metrics.recall_score(y, pr_and_pr[:, 0],
                                      average=None)  # returns recall for each class (3 values, as -1 is class too)
            print "threshold", threshold, "prec", pr[1], "recall", re[-1]
            pp.append([pr[1], re[-1]])  # we take the prec and rec for only virus class

        pp = np.array(pp)
        plt.plot(pp[:, 1], pp[:, 0])
        plt.xlabel("recall")
        plt.ylabel("precision")
        plt.savefig("precision_recall.png")
        plt.clf()
        plt.plot(np.arange(0.5, 0.991, 0.01), pp[:, 0])
        plt.plot(np.arange(0.5, 0.991, 0.01), pp[:, 1])
        plt.xlabel("threshold value")
        plt.ylabel("precision and recall")
        plt.savefig("threshold_precision_recall.png")

    # DROPOUT - notice that in output th
    def dropout_uncertainty(self, X_test, y_test, load_path):
        self.model.load_weights(load_path)
        get_drop_output = K.function([self.model.layers[0].input, K.learning_phase()],
                                     [self.model.layers[-1].output])  # need to define a Keras function

        raw_pred = np.argmax(self.model.predict(X_test, args.batch_size), axis=1)
        truth = y_test
        print metrics.classification_report(truth, raw_pred, target_names=["not virus", "virus"])

        summed_output = np.zeros((len(truth), 2))
        count = np.zeros(len(truth))
        for i in range(1000):
            dr_output = get_drop_output([X_test, 1])[0]
            summed_output += dr_output
            count += np.argmax(dr_output, axis=1)

        evid = np.abs(summed_output[:, 0] - summed_output[:, 1])
        print np.shape(evid), np.max(evid), np.min(evid), np.median(evid)

        # filtered_pred = raw_pred[np.where(evid > 950)[0]
        # filtered_truth = truth[np.where(evid > 950)[0]
        not_sure = np.where(evid < 900)[0]
        pred_th = raw_pred.copy()
        pred_th[not_sure] = 2
        pr = metrics.precision_score(np.concatenate((truth, [2])), np.concatenate((pred_th, [2])), average=None)[1]
        re = metrics.recall_score(np.concatenate((truth, [2])), np.concatenate((pred_th, [2])), average=None)[1]
        # print "filtered_by_evid shape", filtered_pred.shape
        # pr = metrics.precision_score(filtered_truth, filtered_pred, average=None) #returns precision for each class (2 values)
        print "evid ", pr, re
        # print metrics.classification_report(filtered_truth, filtered_pred, target_names=["not virus", "virus","unclassified"])
        print metrics.classification_report(truth, pred_th, target_names=["not virus", "virus", "not_classified"])

        not_sure = np.where(np.abs(count - 500) < 490)[0]
        pred_th = raw_pred.copy()
        pred_th[not_sure] = 2
        # filtered_pred_c = np.concatenate([raw_pred[np.where(count>990)[0]], raw_pred[np.where(count<10)[0]]])
        # filtered_truth_c = np.concatenate((truth[np.where(count>990)[0]], truth[np.where(count<10)[0]]))
        # print "filtered_by_count shape", filtered_pred_c.shape
        # pr=metrics.precision_score(filtered_truth_c, filtered_pred_c, average=None) #returns precision for each class (2 values)
        print "count ", pr, re
        pr = metrics.precision_score(np.concatenate((truth, [2])), np.concatenate((pred_th, [2])), average=None)[1]
        re = metrics.recall_score(np.concatenate((truth, [2])), np.concatenate((pred_th, [2])), average=None)[1]
        print metrics.classification_report(truth, pred_th, target_names=["not virus", "virus", "not_classified"])

    # print metrics.classification_report(filtered_truth_c, filtered_pred_c, target_names=["not virus", "virus","unclassified"])


    def get_weights(self):
        return self.model.get_weights()

    def set_weights(self, weights):
        return self.model.set_weights(weights)


def str2bool(v):  # a helper function that helps read binary parameter values from command line
    return v.lower() in ("yes", "true", "t", "1")


# this part of the code is called only when you run "python network.py ....".
# If you use an object of MetaFF class in another file this will not be called

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--verbose", type=int, choices=[0, 1, 2], default=2)

    parser.add_argument("--layers", type=int, choices=[1, 2, 3, 4, 5], default=2)
    parser.add_argument("--hidden_nodes", type=int, default=1024)
    parser.add_argument("--activation", choices=['tanh', 'relu'], default='relu')
    parser.add_argument("--dropout", type=float, default=0.25)
    parser.add_argument("--batch_size", type=int, default=10)
    parser.add_argument("--epochs", type=int, default=50)
    parser.add_argument("--train_shuffle", choices=['batch', 'true', 'false'], default='true')
    parser.add_argument("--optimizer", choices=['adam', 'rmsprop'], default='adam')
    parser.add_argument("--lr", type=float, default=0.001)
    parser.add_argument("--lr_epochs", type=int, default=1)  # multiply LR with lr_factor after every lr_epochs
    parser.add_argument("--lr_factor", type=float, default=0.95)

    # the next 2 lines are different optinds to deal with unbalanced samples. Use only one.
    parser.add_argument("--class_weight_power", type=float, default=None)
    parser.add_argument("--oversample", type=float,
                        default=None)  # how many times we oversample the less populated class

    # True would mean we use validation set to decide when to stop training
    parser.add_argument("--save_best_model_only", type=str2bool, default=False)

    # we have loads of different datasets, so you need to pick one
    parser.add_argument("--datatype",
                        choices=["avg", "concat", "forward", "clean_forward", "forward_2ORF", "clustered", "gene",
                                 "gene+", "gene-", "gene+-", "newdir"], default="concat")
    parser.add_argument("--binary", type=str2bool, default=True)  # is it binary classification, or 5 classes

    parser.add_argument("--leave_out", type=int,
                        default=None)  # allows to leave out expetiments, based on 2ORF_projects_id_number.txt
    parser.add_argument("--plot_only", type=str2bool,
                        default=False)  # use True for plotting pred,recall,F1 for already existing trained models

    #TODO: input files has to ne requirement if --predict_unknown is true
    parser.add_argument("--predict_unknown", type=str2bool, default=False)
    parser.add_argument("--unknown_rcsu_file", type=file, nargs='?')
    parser.add_argument("--save_predicted_unknown", type=str, nargs='?')
    #parser.add_argument("-f", "--file", dest="filename",help="write report to FILE", metavar="FILE")
    
    parser.add_argument("--train_all_known", type=str2bool, default=False)
    parser.add_argument("--known_rcsu_file", type=file, nargs='?')
    parser.add_argument("--known_value_file", type=file, nargs='?')
    parser.add_argument("--save_path", type=str)  # the most important arguments have no default value. program will crash if no value is given



    args = parser.parse_args()
    print args._get_kwargs

    # args.save_path =  args.save_path+"{epoch:02d}"
    print "using datatype:", args.datatype
    if args.datatype == "avg":
        data = np.loadtxt("average_matrix.txt")
        values = np.loadtxt("average_values.txt", dtype="int")
    elif args.datatype == "gene":
        data = np.loadtxt("data/2genes_RCSU_matrix.txt")
        values = np.loadtxt("data/2genes_values.txt", dtype="int")
    elif args.datatype == "gene+":
        data = np.loadtxt("data/RCSU_plus_matrix.txt")
        values = np.loadtxt("data/RCSU_plus_values.txt", dtype="int")
    elif args.datatype == "gene-":
        data = np.loadtxt("data/RCSU_minus_matrix.txt")
        values = np.loadtxt("data/RCSU_minus_values.txt", dtype="int")
    elif args.datatype == "gene+-":
        data = np.loadtxt("data/gene2+-_matrix.txt")
        values = np.loadtxt("data/gene2+-_values.txt", dtype="int")
    elif args.datatype == "concat":  # DEFAULT#
        data = np.loadtxt("data/after_clust_concat_al2_matrix.txt")
        values = np.loadtxt("data/after_clust_concat_al2_values.txt", dtype="int")
    elif args.datatype == "forward":
        data = np.loadtxt("RCSU59.txt")
        values = np.loadtxt("values.txt", dtype="int")
    elif args.datatype == "clean_forward":
        data = np.loadtxt("data/clean_forward_matrix.txt")
        values = np.loadtxt("data/clean_forward_values.txt", dtype="int")
    elif args.datatype == "forward_2ORF":
        data = np.loadtxt("data/forward_2ORF_matrix.txt")
        values = np.loadtxt("data/forward_2ORF_values.txt", dtype="int")
    elif args.datatype == "clustered":
        M = "data/clustered_rcsu_matrix.txt"
        V = "data/clustered_values_binary.txt"
        print "data from", M, "   and  ", V
        data = np.loadtxt(M)
        values = np.loadtxt(V, dtype="int")
    elif args.datatype == "newdir":
        # TODO: Give them separatelly also
        if args.predict_unknown & args.train_all_known:
            data = np.loadtxt(args.known_rcsu_file)
            values = np.loadtxt(args.known_value_file, dtype="int")
            unknown_rcsu = np.loadtxt(args.unknown_rcsu_file)

    else:  # return error if dataset name is not good
        assert (False)

    nb_options = 2  # number of classes

    if not args.binary:  # special case of 5-classes is treated here
        print "------------ 5 CLASS PREDICTION--------------"
        data = np.loadtxt("RCSU59.txt")
        values = np.loadtxt("values.txt", dtype="int")
        nb_options = 5

    # when we leave projects out we need to cut away the data and leave this as validation set
    if (args.leave_out is not None) and (args.leave_out > -1):
        print "------------LEAVE OUT--------------"
        experiment_nr = np.loadtxt("2ORF_projects_id_number.txt")
        if args.leave_out > np.max(experiment_nr):  # if exp nr out of bounds!!
            args.leave_out = None
            assert (False)
        lo_train_data = data[np.where(experiment_nr != args.leave_out)]
        lo_train_values = values[np.where(experiment_nr != args.leave_out)]
        lo_test_data = data[np.where(experiment_nr == args.leave_out)]
        lo_test_values = values[np.where(experiment_nr == args.leave_out)]
        # usually we just split train and test randomly, but with leave-out we use the separation from above

    if args.leave_out is None:
        X_train, X_test, y_train, y_test = train_test_split(data, values, test_size=0.1, random_state=42)
        ##Save split datasets
        with open("X_train.txt", "w") as self_X_train, \
                open("X_test.txt", "w") as self_X_test, \
                open("Y_train.txt", "w") as self_Y_train, \
                open("Y_test.txt", "w") as self_Y_test:
            np.savetxt(self_X_train, X_train)
            np.savetxt(self_X_test, X_test)
            np.savetxt(self_Y_train, y_train)
            np.savetxt(self_Y_test, y_test)
    else:
        X_train = lo_train_data
        X_test = lo_test_data
        y_train = lo_train_values
        y_test = lo_test_values
        print "we are using EXPERIMENT LEAVE OUT (exp=", args.leave_out, "). \n Train and test sizes:", y_train.shape, y_test.shape

        # print out statistics:
    counts = np.bincount(y_train)
    print "counts of classes in training set :", counts
    counts2 = np.bincount(y_test)
    print "counts of classes in test:", counts2, "ss", np.shape(np.where(y_test == 1))

    # now dealing with oversampling - works only for binary
    if args.oversample is not None:
        print "------------ PERFORMING OVERSAMPLING --------------"
        min_samples = np.min(counts)
        max_samples = np.max(counts)
        smaller_label = np.argmin(counts)
        bigger_label = np.argmax(counts)
        print "Oversampling - originally: smaller_label, min_samples", smaller_label, min_samples

        # we create a new matrix where we will add all the larger class and then the repetitions of smaller
        oversampled_train_X = X_train[np.where(y_train == bigger_label)[0], :]
        oversampled_train_Y = y_train[np.where(y_train == bigger_label)[0]]

        nr_sam = int(args.oversample * min_samples)  # this is what we want in total
        n_repeat = np.repeat(np.where(y_train == smaller_label)[0],
                             int(args.oversample))  # full repetitions of smaller class
        #TODO:
        extra = np.random.choice(np.where(y_train == label)[0], nr_sam % counts[smaller_label])  # randomly picked extra samples
        smaller_samples = np.concatenate((n_repeat, extra))  # (this contains the indexes)
        oversampled_train_X = np.concatenate((oversampled_train_X, X_train[smaller_samples, :]), axis=0)
        oversampled_train_Y = np.concatenate((oversampled_train_Y, y_train[smaller_samples]), axis=0)

        # now we replace the training set with the oversampled training set
        X_train = oversampled_train_X
        y_train = oversampled_train_Y

        print "After oversampling the shapes are", X_train.shape, y_train.shape
        counts = np.bincount(y_train)
        print "counts of classes after oversample:", counts
        counts2 = np.bincount(y_test)
        print "counts of classes in test after oversample (should stay the same):", counts2

    Y_train = np_utils.to_categorical(y_train, nb_options)  # turns it into one-hot vectors
    Y_test = np_utils.to_categorical(y_test, nb_options)  # turns it into one-hot vectors
    assert (np.all(y_train == np.argmax(Y_train, axis=1)))

    # if we use class weights, we want the weight to depend on
    # 1) nr of samples  and 2) a factor that allows us to make the weights more of less penalizing for bigger class
    # in total we multiply the loss with (count_of_this_class/max_count)**power.
    # for bigger class its always 1. for smaller class it is >1, but can be bigger or smaller than (count_of_small/count_of_big)
    if args.class_weight_power is not None:
        print "------------ CLASS WEIGHT CALCULATION --------------"
        counts = np.bincount(np.argmax(Y_train, axis=1)) * 1.0  # *1.0 to make it a float
        args.class_weight_power = (np.max(counts) / counts) ** args.class_weight_power
        print "counts in classw ", args.class_weight_power  # from now on, cl_w_power contains the weights

    # Creating the model - an object of class MetaFF
    model = MetaFF(**vars(args))
    model.init(X_train.shape[1], Y_train.shape[1])

    ############################
    # TODO: remove this later
    print "init vars: ", vars(args)
    print "X_train.shape: %s, Y_train.shape: %s " % (X_train.shape[1], Y_train.shape[1])
    ############################

    # in case we only want to use a model, not train one
    if args.plot_only:
        model.plot_threshold(X_test, y_test, args.save_path + ".hdf5")
        model.dropout_uncertainty(X_test, y_test, args.save_path + ".hdf5")


    # when training
    else:
        model.fit(X_train, Y_train, X_test, Y_test, args.save_path + '.hdf5')
        print "----------------------train set---------------------------"
        report = model.eval(X_train, y_train, args.save_path + ".hdf5")
        print report
        print "----------------------test set---------------------------"
        report = model.eval(X_test, y_test, args.save_path + ".hdf5")
        print report
        model.plot_threshold(X_test, y_test, args.save_path + ".hdf5")
        if args.predict_unknown:
            model.predict_and_save(unknown_rcsu, args.save_path + ".hdf5", args.save_predicted_unknown)