# Function program [ind,k,u] = search(korig,maara,ind,raja)
# ------------------------------------------------------------------
# Search resonances by streching the geometry.
#
# ------------------------------------------------------------------
# CALLS TO : eee.m
# 12/04/00 : Seppo J‰rvemp‰‰ - RNI
# ------------------------------------------------------------------

import numpy as np
    
def search(korig = None,maara = None,ind = None,raja = None): 
    offset1 = - 0.9
    offset2 = 0
    offset3 = 0.9
    laskuri = 1
    k1 = eee(offset1,maara,0)
    k,u,siirto = eee(offset2,maara,0)
    k3 = eee(offset3,maara,0)
    if (ind == 0):
        ero2,ind = np.amin(np.abs(k - korig))
        if (np.abs(k(ind)) < 1e-06):
            ind = ind + 1
        ind
    
    ero1 = korig - k1(ind)
    s1 = np.sign(ero1)
    ero2 = korig - k(ind)
    s2 = np.sign(ero2)
    ero3 = korig - k3(ind)
    s3 = np.sign(ero3)
    #disp([s1 s2 s3])
    if (np.logical_and(np.logical_and((s1 == 1),(s2 == 1)),(s3 == 1))):
        print('?????????????????????????????')
        print(np.array([k1(ind),k(ind),k3(ind)]))
        raise Exception('Can't compress enough')
    
    if (np.logical_and(np.logical_and((s1 == - 1),(s2 == - 1)),(s3 == - 1))):
        jatka = 1
        while (np.logical_and((jatka > 0),(jatka < 4))):

            print('problem, have not stretched enough, trying to fix')
            offset1 = offset1 - 0.9
            k3 = k
            offset3 = offset2
            s3 = s2
            k = k1
            offset2 = offset1
            s2 = s1
            k1 = eee(offset1,maara,0)
            ero1 = korig - k1(ind)
            s1 = np.sign(ero3)
            if (s3 == 1):
                jatka = 0
                print('succeeded')
            else:
                jatka = jatka + 1

        if (jatka != 0):
            raise Exception('FAILED')
    
    while (np.abs(ero2) > raja):

        if (s1 != s2):
            s3 = s2
            offset3 = offset2
            offset2 = 0.5 * (offset1 + offset3)
            s3 = s2
        else:
            if (s2 != s3):
                s1 = s2
                offset1 = offset2
                offset2 = 0.5 * (offset2 + offset3)
                s1 = s2
        k,u,siirto = eee(offset2,maara,0)
        ero2 = korig - k(ind)
        s2 = np.sign(ero2)
        #  disp([k(ind), ero2, offset2, siirto])
        print(np.array([k(ind),ero2,siirto]))

    
    k = k(ind)
    u = u(:,ind)
    # -----------------------------------------------------------------------
    return ind,k,u