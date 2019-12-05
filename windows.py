#!/usr/bin/env python3

# an close_open integer window class, e.g. [0,3) => 0,1,2
class Window:
    def __init__(self, a=0, b=0):
        assert a<=b, f'a={a} <= b={b} does not hold.'
        self.a = a
        self.b = b

    def _empty(self):
        return self.a == self.b
    
    def __add__(self, other):
        if self._empty():
            return other
        if other._empty():
            return self
        new_a = min(self.a, other.a)
        new_b = max(self.b, other.b)
        self.a, self.b = new_a, new_b
        return self
    
    # might return a window or two windows in a tuple
    def __sub__(self, other):
        if not self==other: # not overlap
            return WindowSet(self)
        elif self<=other:
            return WindowSet(self)
        elif self>=other:
            return self._minus(other)
        elif self.a < other.a:
            return WindowSet(Window(self.a, other.a))
        else:
            return WindowSet(Window(other.b, self.b))

    def _minus(self,other):
        assert self>=other and not self.equal(other)
        if self.a==other.a:
            return WindowSet(Window(other.b, self.b))
        elif self.b==other.b:
            return WindowSet(Window(self.a, other.a))
        else:
            return WindowSet(Window(self.a, other.a), Window(other.b, self.b))

    # in operator
    def __contains__(self, val):
        return val>=self.a and val<self.b

    # two windows are equal
    def equal(self, other):
        return self.a==other.a and self.b==other.b

    # two window overlap
    def __eq__(self, other):
        if self._empty() or other._empty():
            return False
        return not (self<other or self>other)

    # overlap, "self" and "other" on two side, overlap in the middle
    def __ne__(self, other):
        flag1 = other.a in self and self.b in other
        flag2 = other.b in self and self.a in other
        return flag1 or flag2

    # not overlaping, "self" on left side of "other"
    def __lt__(self, other):
        return self.b<=other.a
    
    # not overlaping, "self" on right side of "other"
    def __gt__(self, other):
        return self.a>=other.b

    # overlap, "self" is a subset of "other"
    def __le__(self, other):
        return self.a >= other.a and self.b<=other.b

    # overlap, "self" is a superset of "other"
    def __ge__(self, other):
        return self.a<=other.a and self.b>=other.b

    # intersection
    def __and__(self, other):
        if not self==other: # not overlap
            return Window()
        elif self<=other:
            return self
        elif self>=other:
            return other
        elif self.a < other.a:
            return Window(other.a, self.b)
        else:
            return Window(self.a, other.b)


    def __str__(self):
        return f'[{self.a},{self.b})'
    def __repr__(self):
        return f'[{self.a},{self.b})'

# a set of ordered windows
class WindowSet:
    def __init__(self, *win_list):
        self.winset = []
        if not win_list:
            return
        for w in win_list:
            self.append(w,merge=False)
        self.sort()
        self.merge()
    
    # sort the winset list
    def sort(self):
        self.winset.sort(key=(lambda x: x.a))

    # the starting site of windows in list should be in ascending order
    def append(self, w, merge=True):
        if self.winset:
            assert self.winset[-1].a <= w.a
        else:
            self.winset.append(w)
            return

        if not merge:
            self.winset.append(w)
            return

        # merge window "w" into the last element in list 
        curr_win = self.winset[-1]
        if curr_win<w and not curr_win.b==w.a:
            self.winset.append(w)
        elif curr_win.b==w.a:  # in case two window are consecutive
            curr_win.b = w.b
        else:
            curr_win = curr_win + w

    # insert a window into list
    def insert(self,w):
        if not self.winset:
            self.winset.append(w)
            return
        if self.winset[0]>w:
            self.winset[0:0] = [w]
            return
        if self.winset[-1]<w:
            self.winset.append(w)
            return
        for i in range(len(self.winset)-1):
            if self.winset[i]<w and w<self.winset[i+1]:
                self.winset[i+1:i+1] = [w]
                return

    # merge windows into non-overlaps windows
    def merge(self):
        if not self.winset:
            return
        tmpwinset = []
        curr_win = self.winset[0]
        for w in self.winset:
            if curr_win<w and not curr_win.b==w.a:
                tmpwinset.append(curr_win)
                curr_win = w
            elif curr_win.b==w.a:
                curr_win.b = w.b
            else:
                curr_win = curr_win + w
        tmpwinset.append(curr_win)
        self.winset = tmpwinset

    def __str__(self):
        return str(self.winset)


def winset_intersect(winset1: WindowSet, winset2: WindowSet):
    ret_winset = []
    label_vec = []

    if not winset1.winset:
        return winset2.winset,[2]*len(winset2.winset)
    if not winset2.winset:
        return winset1.winset,[1]*len(winset1.winset)

    n1 = len(winset1.winset)
    n2 = len(winset2.winset)
    i1 = 0
    i2 = 0
    w1 = winset1.winset[i1]
    w2 = winset2.winset[i2]
    while i1<n1 or i2<n2:
        if w1<w2:
            ret_winset.append(w1)
            label_vec.append(1)
            i1 += 1
            if i1>=n1:
                w1 = Window(10**6,10**6)
            else:
                w1 = winset1.winset[i1]
        elif w1>w2:
            ret_winset.append(w2)
            label_vec.append(2)
            i2 += 1
            if i2>=n2:
                w2 = Window(10**6,10**6)
            else:
                w2 = winset2.winset[i2]
        elif w1.equal(w2):
            ret_winset.append(w1)
            label_vec.append(3)
            i1 += 1
            i2 += 1
            if i1>=n1:
                w1 = Window(10**6,10**6)
            else:
                w1 = winset1.winset[i1]
            if i2>=n2:
                w2 = Window(10**6,10**6)
            else:
                w2 = winset2.winset[i2]
        elif w1<=w2:
            tmpwinset = w2-w1
            tmplist = [(w1,3)]+[(w,2) for w in tmpwinset.winset]
            tmplist.sort(key=lambda x: x[0].a)
            for x in tmplist[:len(tmplist)-1]:
                ret_winset.append(x[0])
                label_vec.append(x[1])
            i1 += 1
            if i1>=n1:
                w1 = Window(10**6,10**6)
            else:
                w1 = winset1.winset[i1]
            if tmplist[-1][1]==3:
                ret_winset.append(tmplist[-1][0])
                label_vec.append(tmplist[-1][1])
                i2 += 1 
                if i2>=n2:
                    w2 = Window(10**6,10**6)
                else:
                    w2 = winset2.winset[i2]
            else:
                w2 = tmplist[-1][0]
        elif w1>=w2:
            tmpwinset = w1-w2
            tmplist = [(w2,3)]+[(w,1) for w in tmpwinset.winset]
            tmplist.sort(key=lambda x: x[0].a)
            for x in tmplist[:len(tmplist)-1]:
                ret_winset.append(x[0])
                label_vec.append(x[1])
            i2 += 1
            if i2>=n2:
                w2 = Window(10**6,10**6)
            else:
                w2 = winset2.winset[i2]
            if tmplist[-1][1]==3:
                ret_winset.append(tmplist[-1][0])
                label_vec.append(tmplist[-1][1])
                i1 += 1 
                if i1>=n1:
                    w1 = Window(10**6,10**6)
                else:
                    w1 = winset1.winset[i1]
            else:
                w1 = tmplist[-1][0]
        elif w1.a < w2.a:
            tmpw = w1&w2
            ret_winset.append((w1-tmpw).winset[0])
            label_vec.append(1)
            ret_winset.append(tmpw)
            label_vec.append(3)
            w2 = (w2-tmpw).winset[0]
            i1 += 1 
            if i1>=n1:
                w1 = Window(10**6,10**6)
            else:
                w1 = winset1.winset[i1]
        else:
            tmpw = w1&w2
            ret_winset.append((w2-tmpw).winset[0])
            label_vec.append(2)
            ret_winset.append(tmpw)
            label_vec.append(3)
            w1 = (w1-tmpw).winset[0]
            i2 += 1
            if i2>=n2:
                w2 = Window(10**6,10**6)
            else:
                w2 = winset2.winset[i2] 

    return ret_winset, label_vec


def gen_full_win_list(forward_motif_pos_arr, revcom_motif_pos_arr, kmer_len, seq_len):
    """
    Split sequence into consecutive regions of different types
    Input
    forward_motif_pos_arr: start position of forward motif matching on sequence
    revcom_motif_pos_arr: start positio of revcom motif matching on sequence
    kmer_len: length of kmer/motif
    seq_len: length of input sequence
    Output
    win_list: a list of consecutive regions, each element is a tuple in the form of (start pos, end pos)
    win_type_list: a list of region types for corresponding regions in "win_list"
        0 - non motif region
        1 - forward motif region
        2 - revcom motif region
        3 - overlap of forward and revcom region
    """
    forward_winset = WindowSet()
    revcom_winset = WindowSet()
    for pos in forward_motif_pos_arr:
        forward_winset.append(Window(pos, pos+kmer_len))
    for pos in revcom_motif_pos_arr:
        revcom_winset.append(Window(pos,pos+kmer_len))
    reslist, reslabel = winset_intersect(forward_winset, revcom_winset)

    # print()
    # print('forward motif list')
    # print(forward_winset)
    # print('revcom motif list')
    # print(revcom_winset)

    # print()
    # print('intersect win list')
    # for e,n in zip(reslist, reslabel):
    #     print(e, n)

    win_list = []
    win_type_list = []
    curr_pos = 0
    for e,n in zip(reslist, reslabel):
        if curr_pos<e.a:
            win_list.append((curr_pos,e.a))
            win_type_list.append(0)
        curr_pos = e.b
        win_list.append((e.a, e.b))
        win_type_list.append(n)
    if e.b<seq_len:
        win_list.append((e.b, seq_len))
        win_type_list.append(0)

    # print()
    # for e,n in zip(win_list, win_type_list):
    #     print(e, n)

    return win_list,win_type_list


# if __name__=="__main__":
    w1 = Window(1,3)
    w2 = Window(2,4)
    w3 = Window(5,8)
    w4 = Window(1,8)
    print(w1)
    print(w1!=w2)
    print(w1+w2)
    print(w1!=w3)
    print(w1!=w4)
    print(w1<=w4)
    print(w2==w4)
    print(w4-w2)
    # for w in w4-w2:
    #     print(w)

    w1 = Window(1,3)
    w2 = Window(6,8)
    w3 = Window(10,13)

    w4 = Window(2,4)
    w5 = Window(5,12)
    w6 = Window(14,18)
    w7 = Window(20,28)

    # ws1 = WindowSet(w1,w2,w3)
    ws2 = WindowSet(w4,w5,w6,w7)
    ws1 = WindowSet()

    print()
    print(ws1)
    print(ws2)

    reslist, reslabel = winset_intersect(ws1, ws2)

    print('res')
    # for e,n in zip(reslist, reslabel):
    #     print(e, n)
    
    curr_pos = 0
    for e,n in zip(reslist, reslabel):
        if curr_pos<e.a:
            print(Window(curr_pos,e.a),0)
        curr_pos = e.b
        print(e,n)

    # print(reslist)
    # print(reslabel)

    forward_motif_pos_arr = [3, 8, 17, 20]
    # revcom_motif_pos_arr = [5, 12, 17, 19, 30]
    revcom_motif_pos_arr = []

    kmer_len = 4
    seq_len = 50

    win_list, win_type_list = gen_full_win_list(forward_motif_pos_arr, revcom_motif_pos_arr, kmer_len, seq_len)




