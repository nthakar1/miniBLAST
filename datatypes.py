
# Not actually using these currently but we'll see

class RefSeqNode:
    def __init__(self, kmer):
        self.kmer = kmer
        self.pos = None
        self.left = None
        self.right = None

class BST:
    def __init__(self):
        self.root = None

    def buildTreeFromDictionary(self, d):
        keys = sorted(d.items())
        self.root = self._build(keys, 0, len(d)-1)

    def _build(self, d, leftIndex, rightIndex):
        if leftIndex > rightIndex:
            return None
        
        midIndex = (leftIndex + rightIndex) // 2
        node = RefSeqNode(d[midIndex])

        node.left = self._build(d, leftIndex, midIndex-1)
        node.right = self._build(d, midIndex+1, rightIndex)

        return node