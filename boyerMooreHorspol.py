from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from itertools import groupby
import sys

#------------MAIN METHOD--------------------
#reads in fasta file, runs program to either count or find pos
def main():
  if (checkArgs(sys.argv)):
    faFile = sys.argv[1]
    marker = sys.argv[2]
    flag = sys.argv[3]
    alphaSize = len(''.join(set(marker)))
    for fastaRec in SeqIO.parse(faFile, 'fasta'):
      table = Table(alphaSize, marker)
      preProcessTable(table, marker)
      run(flag, table, marker, fastaRec.seq)


#---------METHODS------------------------

#runs either occurence counter or position finder
def run(flag, table, marker, seq):
  if (int(flag) == 0):
    countOccurences(table, marker, seq)
  elif (int(flag) == 1):
    findPos(table, marker, seq)

#simply checks arguements entered (this could be edited to created more arg checks)
def checkArgs(args):
  if (len(args) == 4):
    if (int(args[3]) == 0 or int(args[3]) == 1):
      if (len(args[2]) >= 2):
        return True
      else:
        print "Marker must be at least 2 characters"
        return False
    else:
      print "Last Arguement must be 0 or 1"
      return False
  else:
    return False

#sets the values in the table
def preProcessTable(table, marker):
  for index in range(table.size):
    table.matrix[index][0] = table.matrix[index][0] - 1
  for index in range(len(marker) - 2):
    couple =  marker[index:index+2]
    x = table.set.index(couple[0])
    y = table.set.index(couple[1])
    table.matrix[x][y] = len(marker) - (index+2)

#counts the amount of times marker occurs in seq
def countOccurences(table, marker, seq):
  occurence = 0
  index = 0
  while index < (len(seq) - len(marker)):
    partition = seq[index:(index+len(marker))]
    if (marker in partition):
      occurence += 1
    couple = seq[index+len(marker)-2:index+len(marker)]
    if (couple[0] in table.set and couple[1] in table.set):
      x = table.set.index(couple[0])
      y = table.set.index(couple[1])
      shift = table.matrix[x][y]
      index = index + shift
    else:
      index = index + len(marker)
  print "The Marker: ", marker, " occurs: ", occurence

#finds the start and end positions of the marker
def findPos(table, marker, seq):
  index = 0
  while index < (len(seq) - len(marker)):
    partition = seq[index:(index+len(marker))]
    if (marker in partition):
      print "Start position: ", index+1, " End Position: ", index+len(marker)+1
    couple = seq[index+len(marker)-2:index+len(marker)]
    if (couple[0] in table.set and couple[1] in table.set):
      x = table.set.index(couple[0])
      y = table.set.index(couple[1])
      shift = table.matrix[x][y]
      index = index + shift
    else:
      index = index + len(marker)

#-----CLASSES---------------------CONTAINS: Table
#uses an alphabet to create a table with shift count
class Table(object):
  def __init__(self, size=None, marker=None):
    self.matrix = None
    self.size = size
    self.set = ''.join(k for k, g in groupby(sorted(marker)))
    length = len(marker)
    if (size != None and size > 0):
      self.matrix = [[length for x in range(size)] for x in range(size)] 

  def printTable(self):
    for col in range(self.size):
      print self.matrix[col][0:self.size]

main()
