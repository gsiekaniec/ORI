#!/usr/bin/env python3
# coding: utf-8

def vraissemblance(strains,nb_strains,reads,abondances,length):
    '''Computation of the likelihood'''
    v = 0
    res = 0
    abondance_max = 0
    for i in range(nb_strains):
        abondance_max = abondance_max+(abondances[i]*length[i])
    for strs in reads:
        for num in range(nb_strains):
            cij = 0
            if strains[num] in strs:
                cij = 1
            actualstrain = (abondances[num]*length[num])/abondance_max
            res = res+(actualstrain*cij)
        if v == 0:
            v = res
        else:
            v = v*res
    return v

def expectation(strains,nb_strains,reads,abondances):
    '''Computation of the expectation part of the EM algorithm'''
    res = []
    for num in range(nb_strains):
        nj = 0
        for strs in reads:
            abondance_max = 0
            for i in range(nb_strains):
                if strains[i] in strs:
                    abondance_max = abondance_max+abondances[i]
            if strains[num] in strs:
                nj = nj+ (abondances[num]/abondance_max)
        res.append(nj)
    return res        
     
def maximization(strains,nb_strains,reads,njs,length):
    '''Computation of the maximization part of the EM algorithm'''
    res = []
    nk_max = 0
    for i in range(nb_strains):
        nk_max = nk_max+(njs[i]/length[i])
    for num in range(nb_strains):
        aj = (njs[num]/length[num])/nk_max
        res.append(aj)
    return res

