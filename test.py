import numpy as np
from BasisState import BasisState
from DisorderParameter import KitaevDisorderParameter
from KitaevHamiltonianBlock import NaiveKitaevHamiltonianBlock

# naive_computed_element
N = 6
Q = 4

J = KitaevDisorderParameter(N, 1., False)
H = NaiveKitaevHamiltonianBlock(N, Q, J)

ar0123 = [0,1,2,3]
ar0145 = [0,1,4,5]

bra = BasisState(_indices=ar0145, _coefficient=4)
ket = BasisState(_indices=ar0123, _coefficient=4)

left = H.matrix[bra.get_state_number(), ket.get_state_number()]
right_answer = -np.sqrt(2.) / pow(N, 3./2.) * J.elem(4,5,2,3)
assert np.isclose(left, right_answer), (left, right_answer)

