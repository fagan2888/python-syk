from typing import Union, List

import numpy as np
from scipy.special import binom as binomial
#include <sstream>
#include <iostream>
#include <assert.h>
#include "BasisState.h"
#include "FockSpaceUtils.h"

# Constructs the object with the vacuum state
class BasisState:
    def __init__(
            self, _indices: Union[List[int], int] = None,
            _num_indices: int = None,
            _coefficient: int = None,
            state_number: int = None, Q: int = None,
            global_state_number: int = None,
    ):
        if _indices is not None and _coefficient is not None:
            self.indices = _indices
            self.coefficient = _coefficient
        #elif _indices is not None and _num_indices is not None and _coefficient is not None:
        #    self.indices = list<int>(_indices, _indices + _num_indices)
        #    self.coefficient = _coefficient
        elif state_number is not None and Q is not None:
            self.indices = state_number_to_occupations(state_number, Q)
            self.coefficient = 1
        elif global_state_number is not None:
            self.indices = []
            assert global_state_number >= 0
            i = 0
            while global_state_number != 0:
                bit = global_state_number % 2
                global_state_number /= 2
                if bit:
                    self.indices.append(i)
                i += 1
            self.coefficient = 1
        else:
            self.coefficient = 1

    # Act with c_i
    def annihilate(self, i: int) -> None:
        if self.coefficient == 0:
            return
        index_found = False

        for idx, _iter in enumerate(self.indices):
            if _iter == i:
                # Found the position
                index_found = True
                break
            elif _iter < i:
                # We need to drag c_i across this, so we flip the sign
                self.coefficient *= -1
            else:
                # No need to go any further, index not found in state
                # The index i doesn't appear in the state, so c_i kills it
                break
        if index_found:
            del self.indices[idx]
        else:
            self.coefficient = 0

    # Act with c^\dagger_i
    def create(self, i: int) -> None:
        if self.coefficient == 0:
            return
        for idx, _iter in enumerate(self.indices):
            if _iter == i:
                # Acting twice with same creation operator kills the state
                self.coefficient = 0
                return
            elif _iter < i:
                # We need to drag c^\dagger_i across this, so we flip the sign
                self.coefficient *= -1
            else:
                # We found the position
                break

        self.indices.insert(idx, i)

    def get_state_number(self) -> int:
        return occupations_to_state_number(self.indices)

    def get_global_state_number(self) -> int:
        return sum(pow(2, i) for i in self.indices)

    def charge(self) -> int:
        assert self.coefficient != 0
        return len(self.indices)

    def is_zero(self) -> bool:
        return self.coefficient == 0

    def next_state(self) -> 'BasisState':
        return BasisState(
            state_number=self.get_state_number() + 1,
            Q=self.charge()
        )

#BasisState& BasisState::operator = (const BasisState& other) {
#    if (this != &other) {
#        indices = other.indices
#        coefficient = other.coefficient
#    }
#
#    return *this
#}

    def __eq__(self, other: 'BasisState') -> bool:
        return (self.indices == other.indices) and (self.coefficient == other.coefficient)

    def __ne__(self, other: 'BasisState') -> bool:
        return not (self == other)

    def __str__(self):
        if self.coefficient == 0:
            return "0"
        else:
            s = ""
            if self.coefficient == -1:
                s += "-"
            elif self.coefficient != 1:
                s += str(self.coefficient) + "*"
            s += "|"
            for _iter in self.indices:
                s += str(_iter) + ","
            s = s[:-1]  # remove the last ","
            s += ">"
            return s

# Find maximal c_k such that (c_k choose k) <= N
def find_maximal_ck(state_number: int, k: int) -> int:
    ck = 0
    while int(binomial(ck, k)) <= state_number:
        ck += 1
    return ck - 1

def state_number_to_occupations(state_number: int, Q: int) -> List[int]:
    k = Q
    indices = []
    while k > 0:
        ck = find_maximal_ck(state_number, k)
        indices.insert(0, ck)
        state_number -= int(binomial(ck, k))
        k -= 1
    assert state_number == 0
    return indices

def occupations_to_state_number(indices: List[int]) -> int:
    k = 1
    result = 0
    for _iter in indices:
        result += int(binomial(_iter, k))
        k += 1
    return result

#GlobalStateIterator::GlobalStateIterator(const Space& _space) {
#    space = _space
#    global_state = 0
#    parity_ordered_state = 0
#
#    for (int Q = 0 Q <= space.Nd Q++) {
#        seen_states_by_charge[Q] = 0
#    }
#}
#
#bool GlobalStateIterator::done() {
#    return global_state >= space.D
#}
#
#void GlobalStateIterator::next() {
#    if (done()) return
#    global_state++
#    if (done()) return
#
#    int state_Q = __builtin_popcount(global_state)
#
#    # Offset within charge sector
#    parity_ordered_state = seen_states_by_charge.at(state_Q)
#
#    if (state_Q % 2 == 0) {
#        # Add all previous even sectors
#        for (int Q = 0 Q < state_Q Q += 2) {
#            parity_ordered_state += Q_sector_dim(space.Nd, Q)
#        }
#    }
#    else {
#        # For odd charges, add all even sectors, and all previous
#        # odd sectors
#        for (int Q = 0 Q <= space.Nd Q += 2) {
#            parity_ordered_state += Q_sector_dim(space.Nd, Q)
#        }
#
#        for (int Q = 1 Q < state_Q Q += 2) {
#            parity_ordered_state += Q_sector_dim(space.Nd, Q)
#        }
#    }
#
#    assert(parity_ordered_state >= 0)
#    assert(parity_ordered_state < space.D)
#
#    seen_states_by_charge.at(state_Q)++
#}
