#------------ALISSON BARIL--------------------

import Bio
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import Phylo
from itertools import groupby
import sys

#------------VARIABLES--------------------
cost = [[0, 2, 1, 2],
        [2, 0, 2, 1],
        [1, 2, 0, 2],
        [2, 1, 2, 0]]

dictionary = {}

seqs = {};

taken = []

#------------MAIN METHOD--------------------
#reads in fasta file and the newick file
def main():
  if (checkArgs(sys.argv)):
    faFile = sys.argv[2]
    tree = sys.argv[1]
    scores = Root()
    size = 0
    for fastaRec in SeqIO.parse(faFile, 'fasta'):
      seqs[fastaRec.id] = fastaRec.seq
      if (size is 0):
        size = len(fastaRec.seq)
    for i in range(size):
      scores = run(tree, i, scores)
    scores.printRoots()


#---------METHODS------------------------

#simply checks arguements entered (this could be edited to created more arg checks)
def checkArgs(args):
  if (len(args) == 3):
    return True
  else:
    return False

#method creates terminal score tables and prepare to calc all tables
def run(tree, index, scores):
  tree = Phylo.read(tree, "newick")
  taken = tree.get_terminals()
  for leaf in taken:
    dictionary[leaf.name] = newLeafScore(seqs[leaf.name], index)
  nonterm = []
  for node in tree.get_nonterminals(order='level'):
    nonterm.append(node.name)
  loopThrough(nonterm, tree, taken)
  lastNode = taken.pop()
  scores.addRootScore(dictionary[lastNode.name])
  dictionary.clear()
  taken[:] = []
  return scores

#recursively loops through calcing score tables
def loopThrough(nodes, tree, taken):
  if (len(nodes) > 0):
    cladeList = []
    for i in range(0, len(dictionary)):
      for j in range(1, len(dictionary)):
        if (i is not j):
          cladeList = tree.trace(taken[i], taken[j])
          if (len(cladeList) == 2 and cladeList[0].name in nodes):
            dictionary[cladeList[0].name] = newScore(taken[i], taken[j])
            taken.append(cladeList[0])
            nodes.remove(cladeList[0].name)
    loopThrough(nodes, tree, taken)
  else:
    return        
 

#creates score tables for the leafs
def newLeafScore(seq, index):
  a, c, g, t = float("inf"), float("inf"), float("inf"), float("inf")
  if (seq[index] is "A"):
    a = 0
  if (seq[index] is "C"):
    c = 0
  if (seq[index] is "G"):
    g = 0
  if (seq[index] is "T"):
    t = 0
  return Score(a, c, g, t)

#creates score table for node
def newScore(clade1, clade2):
  scorelist = []
  list1, list2 = [], []
  for i in range(4):
    list1[:] = []
    list2[:] = []
    for j in range(4):
      list1.append(cost[i][j] + dictionary[clade1.name].scoreByNum(j))
      list2.append(cost[i][j] + dictionary[clade2.name].scoreByNum(j))
    scorelist.append(min(list1) + min(list2))
  return Score(scorelist[0], scorelist[1], scorelist[2], scorelist[3])
    

#-----CLASSES---------------------CONTAINS: Score, Root
#Score is a simple class that keeps track of the a,c,g,t vals
class Score(object):
  def __init__(self, a = None, c = None, g = None, t = None):
    self.a = a
    self.c = c
    self.g = g
    self.t = t 

  def printScore(self):
    print (self.a, self.c, self.g, self.t)

  def scoreByNum(self, num = 0):
    if (num is 0):
      return self.a
    if (num is 1):
      return self.c
    if (num is 2):
      return self.g
    else:
      return self.t

#an array of Scores made for each Root 
class Root(object):
  def __init__(self):
    self.rootScores = []

  def addRootScore(self, score = None):
    self.rootScores.append(score)

  def printRoots(self):
    aList, cList, gList, tList = "", "", "", ""
    for a in self.rootScores:
      aList = aList + str(a.scoreByNum(0)) + " "
    for c in self.rootScores:
      cList = cList + str(c.scoreByNum(1)) + " "
    for g in self.rootScores:
      gList = gList + str(g.scoreByNum(2)) + " "
    for t in self.rootScores:
      tList = tList + str(t.scoreByNum(3)) + " "
    print("A: ", aList)
    print("C: ", cList)
    print("G: ", gList)
    print("T: ", tList)

main()
