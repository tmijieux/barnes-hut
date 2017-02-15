#ifndef BH_PROC_H
#define BH_PROC_H

#include <iostream>

namespace barnes_hut {

class proc {
private:
    int _rank;
    int _group_size;

public:
    inline int rank() const { return _rank; }
    inline int size() const { return _group_size; }
    void print() const;
    void print(std::ostream&) const;
    proc();
    ~proc();
};


}; // end namespace barnes_hut

#endif // BH_PROC_H
