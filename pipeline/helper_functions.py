import numpy as np

# datetime format - should probably read this from a config file and not hard coded here
_datetimeformat_ymd = '%y%m%d'
_datetimeformat_ydm = '%y%d%m'

def Get1FromNestedArray(nestedArray):
    if nestedArray.size == 0: return None
    unpackedVal = nestedArray
    while unpackedVal.shape != (): 
        if unpackedVal.size == 0: return None
        unpackedVal = unpackedVal[0]
    return unpackedVal

def GetListFromNestedArray(nestedArray):
    if nestedArray.size == 0: return None
    unpackedVal = nestedArray
    l = []
    if unpackedVal.size == 0: return None
    for j in np.arange(unpackedVal.size):
        tmp = unpackedVal.item(j)
        try: tmp = Get1FromNestedArray(tmp)
        except: pass
        l.append(tmp)            
    return l    