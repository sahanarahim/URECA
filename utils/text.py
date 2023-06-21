'''
Contains functions used to generate the text summary of a network.
'''

def make_text(elements):    
    '''Given all edges in the KnowledgeNet, it makes the text summary'''
    pubmedLink = '<span class = "pubmed-link" data-pubmed-id = "%s" data-source = "%s" data-typa = "%s" data-target = "%s">%s</span>'
    #Paragraph order: node by the highest degree. Sentence order in a paragraph: node,interaction type by the number of targets.  
    topicDic = {}
    nodeDegree, nodeSentenceDegree = {}, {}
    for i in elements:
        if i.id not in topicDic:
            topicDic[i.id]={}
            nodeDegree[i.id]=0
            nodeSentenceDegree[i.id]= {}
        if i.inter_type not in topicDic[i.id]:
            topicDic[i.id][i.inter_type] = [[i.target, i.publication]]
            nodeSentenceDegree[i.id][i.inter_type] = 1
            nodeDegree[i.id] += 1
        else:
            topicDic[i.id][i.inter_type] += [[i.target, i.publication]]
            nodeDegree[i.id] += 1
            nodeSentenceDegree[i.id][i.inter_type] += 1

    sorted_nodes = sorted(nodeDegree, key=lambda x: nodeDegree[x], reverse=True)    
    # print(sorted_nodes)
    save = []
    for i in sorted_nodes:
        sorted_sentences = sorted(nodeSentenceDegree[i], key=lambda x: nodeSentenceDegree[i][x], reverse=True)
        tempSentences = []
        for j in sorted_sentences:
            text = '<span  style="color: red;">' + i + '</span>' + ' ' + '<span  style="color: orange;">' + j + '</span>' + ' '
            temp = {}
            tempRefs = []
            for k in topicDic[i][j]:
                if k[0] not in temp:
                    temp[k[0]] = [pubmedLink % (k[1], i, j, k[0], k[1])]
                else:
                    temp[k[0]] += [pubmedLink % (k[1], i, j, k[0], k[1])]
            
            for target in temp:
                tempRefs += [target + ' (' + ', '.join(list(set(temp[target])))+ ')']
            tempSentences.append(text + ', '.join(tempRefs))

        finishedSentence= '. '.join(tempSentences) + '.'
        save.append("<div id = \""+i+"\">" + finishedSentence + "</div>")
    return '<br><br>'.join(save)