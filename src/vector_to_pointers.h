#ifndef VECTOR_TO_POINTERS_H
#define VECTOR_TO_POINTERS_H

#include <vector>

template<typename T, class V>
std::vector<T*> vector_to_pointers(std::vector<V>& input) {
    std::vector<T*> output(input.size());
    auto oIt = output.begin();
    for (auto& i : input) {
        *oIt = static_cast<T*>(i.begin());
        ++oIt;
    }
    return output;
}

template<typename T, class V>
std::vector<const T*> vector_to_pointers(const std::vector<V>& input) {
    std::vector<T*> output(input.size());
    auto oIt = output.begin();
    for (auto& i : input) {
        *oIt = static_cast<const T*>(i.begin());
        ++oIt;
    }
    return output;
} 

#endif
