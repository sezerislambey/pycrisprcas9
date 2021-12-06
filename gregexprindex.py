from Bio.Seq import Seq
def gregexpr_index(pattern, string):

    index = []
    for i in range(len(pattern)):
        code_index = string.find(pattern[i])

        while code_index >= 0:
        
            index.append(code_index)
            code_index = string.find(pattern[i], code_index+1)
            if code_index == -1:
                break

    return index
