import sys

def traverse(invertedindex, checker, elem, clusterelem):
    stk = []
    stk.append(elem)
    while len(stk) != 0:
        elem = stk[len(stk) - 1]
        del stk[len(stk) - 1]
        if checker.has_key(elem) == False:
            checker[elem] = True
            clusterelem.append(elem)
            index = invertedindex[elem]
            for elem2 in index:
                stk.append(elem2)

argvs = sys.argv

f = open(argvs[1])
line = f.readline()

maxid = 0
while line:
    tmp = line.replace('\n', '').split(" ")

    if maxid < int(tmp[0]):
        maxid = int(tmp[0])
    if maxid < int(tmp[1]):
        maxid = int(tmp[1])
    line = f.readline()
f.close

invertedindex = []
for i in (range(maxid + 1)):
    invertedindex.append([])

f = open(argvs[1])
line = f.readline()

while line:
    tmp = line.replace('\n', '').split(" ")
    invertedindex[int(tmp[0])].append(int(tmp[1]))
    invertedindex[int(tmp[1])].append(int(tmp[0]))
    line = f.readline()

#print(maxid)
#print

independent = []
clusters = []
checker = {}
for i in (range(maxid + 1)):
    clusterelem = []
    traverse(invertedindex, checker, i, clusterelem)
#    print(str(i) + " " + str(len(clusterelem)))
    if len(clusterelem) > 1:
        clusters.append(clusterelem)
    elif len(clusterelem) == 1:
        independent.append(clusterelem[0])
    

cid = 0        
for cluster in clusters:
    print(str(cid) + " #" + str(len(cluster)))
    cid += 1
    for c in cluster:
        print(str(c)) , 
    print

print("isolated cluster: #" + str(len(independent)))
for elem in independent:
    print(elem) ,
print

    # int(tmp[0]), int(tmp[1]), float(tmp[2])
    
#f = open(argvs[1])
#line = f.readline()
#while line:
#    line = f.readline()
#    tmp = line.replace('\n', '').split(" ")

    

#    print float(tmp[2])

