from math import ceil, log2

Input = {
  'a': 10,
  'e': 15,
  'i': 12,
  'o': 3,
  'u': 4,
  's': 13,
  't': 1
}

#  Supporting Functions
class BinaryTreeNode:
    def __init__(self, data):
        self.data = data
        self.leftChild = None
        self.rightChild = None

def Entropy(Input_Dict):
    TotalChar = 0
    H = 0
    for key in Input_Dict:
        TotalChar = TotalChar + Input_Dict[key]
    for key in Input_Dict:
        p = Input_Dict[key] / TotalChar
        H = H + p*log2(1/p)
    return H

def AvgCodeLen(Input_Dict,Huff_Code_Dict):
    TotalChar = 0
    L = 0
    for key in Input_Dict:
        TotalChar = TotalChar + Input_Dict[key]
    for key in Input_Dict:
        p = Input_Dict[key] / TotalChar
        CodeLen = len(Huff_Code_Dict[key])
        L = L + p*CodeLen
    return L

def Efficiency(Hm,L):
    if L != 0:
        return (Hm/L)*100
    else:
        return 100

# Data Compression Algorithms
def ShannonFano(Data):
    # Tree based algo
    pass

def HuffmanCode(Input_Dict):
    # Tree based algo
    Huff_Code_Dict = {}
    Input_Dict = sorted(Input_Dict.items(), key=lambda x: x[1])

    for key in Input_Dict:
        pass

def LZW(Data):
    UniqueElementsCount = len(set(Data))
    MaxLen = ceil(log2(UniqueElementsCount))
    HelperSet = set(Data)
    k = 0
    for Element in HelperSet:
        BinNum = str(bin(k).replace("0b", ""))
        Data = [BinNum if Index == Element else Index for Index in Data]

    k = 0
    PartitonedData = []
    TempString = Data[k]
    while True:
        if not(TempString in PartitonedData):
            PartitonedData.append(TempString)
            k = k + 1
            if (k < len(Data)):
                TempString = Data[k]
            else:
                break
        else:
            k = k + 1
            if (k < len(Data)):
                TempString = TempString + Data[k]
            else:
                PartitonedData.append(TempString)
                break
    return PartitonedData

# Error Correcting Algorithms
def GolayCode(Data):
    pass

def HammingCode(Data,Parity):
    m = len(Data)
    r = 0
    Parity = int(Parity)
    while (pow(2,r) < (m + r +1)):
        r = r + 1
    Index = 1
    k = 1
    HamMess = []
    while (Index <= (m+r)):
        if (Index & (Index-1) == 0):
            HamMess.append(2)
        else:
            HamMess.append(int(Data[m-k]))
            k = k + 1
        Index = Index + 1
    HamMess.reverse()
    for Index in range(1,r+1):
        SubIndex = int(pow(2,Index-1))
        ParityCounter = 0
        NumCounter = 0
        while SubIndex <= (m + r):
            if (NumCounter < int(pow(2,Index-1))):
                if HamMess[m + r - SubIndex] == 1:
                    ParityCounter = ParityCounter + 1
                NumCounter = NumCounter + 1
                Instance = SubIndex
            else:
                if (SubIndex == (Instance + NumCounter)):
                    NumCounter = 0
            SubIndex = SubIndex + 1
        if Parity == 0:
            if ParityCounter%2 == 0:
                HamMess[m + r - int(pow(2,Index-1))] = 0
            else:
                HamMess[m + r - int(pow(2,Index-1))] = 1
        else:
            if ParityCounter%2 == 0:
                HamMess[m + r - int(pow(2,Index-1))] = 1
            else:
                HamMess[m + r - int(pow(2,Index-1))] = 0
    HamMess = ''.join(map(str,HamMess))
    return HamMess

def BCHCode(Data):
    pass

def ReedSolomonCode(Data):
    pass

def CRC(Data,CRCGen):
    Padding = ''
    for i in range(1,len(CRCGen)):
        Padding = Padding + '0' 
    ExtendedData = Data + Padding
    Index = 1
    Start = 0
    while True:
        for Index in range(0,len(CRCGen)):
            if (Index + Start) < len(ExtendedData):
                ExtendedData[Index+Start] = int(ExtendedData[Index+Start])^int(CRCGen[Index])
        for Index in range(0,len(CRCGen)):
            Start = Index
            if int(ExtendedData[Index+Start]) == 1:
                break
        if (Start + len(CRCGen) - 1) > len(ExtendedData):
            # CRCRemainder = ??
            break

    return CRCRemainder

# Driver Code
Data = input()
LZWData = LZW(Data)
print(LZWData)

# Data = input()
# CRCGen = input()
# CRCBits = CRC(Data,CRCGen)
# print(CRCBits)

# Data = input()
# Parity = int(input())
# HamCode = HammingCode(Data,Parity)
# print(HamCode)

# Huff_Code = HuffmanCode(Input)
# Hm = Entropy(Input)
# L = AvgCodeLen(Input,Huff_Code)
# Eta = Efficiency(Hm,L)

# print(Huff_Code)
# print(Hm)
# print(L)
# print(Eta)

# node1 = BinaryTreeNode(50)
# node2 = BinaryTreeNode(20)
# node3 = BinaryTreeNode(45)
# node4 = BinaryTreeNode(11)
# node5 = BinaryTreeNode(15)
# node6 = BinaryTreeNode(30)
# node7 = BinaryTreeNode(78)

# node1.leftChild = node2
# node1.rightChild = node3
# node2.leftChild = node4
# node2.rightChild = node5
# node3.leftChild = node6
# node3.rightChild = node7

Raw = input()