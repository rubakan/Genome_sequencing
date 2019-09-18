#!/usr/bin/env python3

#*************************************************************************
#
#   Program:    Titanic Prediction
#   File:       Titanic_predict.py
#   
#   Version:    Python 3
#   Date:       27.06.2019
#   
#   
#   Copyright:  (c) NMBU / Rubakan Thurupan
#   Author:     Rubakan Thurupan
#   Address:    Faculty of Chemistry, Biotechnology and Food Science (KBM),
#               Norwegian University of Life Science,
#               Universitetstunet 3,
#               1430 Ås.
#               Norway
#              
#   EMail:      rubakan.thurupan@nmbu.no

#*************************************************************************************************************************
#
#
#
#
#********************************************************************************************************************
#
#   Description:
#   ============
#
#********************************************************************************************************************
#
#   Usage:
#   ======
#   Extract a Excel file to Python Script.
#
#********************************************************************************************************************
#
#   Revision History:
#   =================
#   V1.0   17.09.19  Original
#   
#                                    
#                    
#
#
#********************************************************************************************************************

import pandas as pd

train = pd.read_csv('/Users/Rubakan/git/machine_learning/feature_engineering/Titanic_Test.csv')
test = pd.read_csv('/Users/Rubakan/git/machine_learning/feature_engineering/Titanic_train.csv')

#Drop features we are not going to use
train = train.drop(['Name','SibSp','Parch', 'Ticket', 'Fare', 'Cabin', 'Embarked'],axis=1)
test = test.drop(['Name','SibSp','Parch', 'Ticket', 'Fare', 'Cabin', 'Embarked'],axis=1)

#Look at the first 3 rows of our training data
train.head(3)



for df in [train,test]:
    df['Sex_binary']=df['Sex'].map({'male':1,'female':0})
    
#Fill in missing age values with 0 (presuming they are a baby if they do not have a listed age)
train['Age'] = train['Age'].fillna(0)
test['Age'] = test['Age'].fillna(0)

#Select feature column names and target variable we are going to use for training
features = ['Pclass','Age','Sex_binary']
target = 'Survived'

#Look at the first 3 rows (we have over 800 total rows) of our training data.; 
#This is input which our classifier will use as an input.
train[features].head(3)

train[target].head(3).values

from sklearn.tree import DecisionTreeClassifier

#Create classifier object with default hyperparameters
clf = DecisionTreeClassifier()  

#Fit our classifier using the training features and the training target values
clf.fit(train[features],train[target]) 
