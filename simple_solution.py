dt = 0.001
epsilon = 0.00001

def equilibrate(k_a, k_d, c0_E, c0_L, c0_EL) -> (float, float, float):
    delta = 1
    c_E = c0_E
    c_L = c0_L
    c_EL = c0_EL
    i = 0
    while delta > epsilon:
        delta_a = dt * k_a * c_E * c_L
        delta_d = dt* k_d * c_EL
        delta = delta_a - delta_d
        c_EL += delta
        c_E -= delta
        c_L -= delta
        
        print(i, c_L, c_E, c_EL)
        i += 1

    return c_L, c_E, c_EL
    
c_L, c_E, c_EL = equilibrate(1.0, 1.0, 1.0, 1.0, 0.0)

print(c_L, c_E, c_EL )
    