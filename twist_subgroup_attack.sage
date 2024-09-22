from sage.all import EllipticCurve, GF, factor, crt, gcd, discrete_log, power_mod, ZZ
from itertools import permutations

# Define the secp256k1 curve in SageMath
p = 2**256 - 2**32 - 977  # secp256k1 prime
n = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141  # Order of the curve
a = 0
b = 7
secp256k1 = EllipticCurve(GF(p), [a, b])

# Use provided public key coordinates
pubx = ZZ("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx", 16)
puby = ZZ("123abcblablablaxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx", 16)
public_key = (pubx, puby)
print(f"Provided Public Key: {public_key}")

# Verify if the provided public key coordinates are on the secp256k1 curve
if secp256k1.is_on_curve(pubx, puby):
    print("The provided public key is valid on the secp256k1 curve.")
else:
    raise ValueError("The provided public key is not valid on the secp256k1 curve.")

# Define the twist curves
curves = {
    "Twist Curve 1": EllipticCurve(GF(p), [0, -7]),
    "Twist Curve 2": EllipticCurve(GF(p), [0, -6]),
    "Twist Curve 3": EllipticCurve(GF(p), [0, -1]),
    "Twist Curve 4": EllipticCurve(GF(p), [0, 17])
}

# Function to check the Legendre symbol
def legendre_symbol(a, p):
    if a == 0:
        return 0
    ls = power_mod(a, (p - 1) // 2, p)
    return -1 if ls == p - 1 else ls

# Function to find square roots modulo a prime using the Tonelli-Shanks algorithm
def tonelli_shanks(n, p):
    assert legendre_symbol(n, p) == 1, "n is not a quadratic residue modulo p"
    if p % 4 == 3:
        return power_mod(n, (p + 1) // 4, p)
    q = p - 1
    s = 0
    while q % 2 == 0:
        q //= 2
        s += 1
    z = 2
    while legendre_symbol(z, p) != -1:
        z += 1
    m = s
    c = power_mod(z, q, p)
    t = power_mod(n, q, p)
    r = power_mod(n, (q + 1) // 2, p)
    while t != 0 and t != 1:
        t2i = t
        for i in range(1, m):
            t2i = (t2i * t2i) % p
            if t2i == 1:
                break
        b = power_mod(c, 1 << (m - i - 1), p)
        m = i
        c = (b * b) % p
        t = (t * c) % p
        r = (r * b) % p
    return r

# Function to transform and validate points on twist curves
def find_valid_y_on_twist(x, curve):
    rhs = x**3 + curve.a4() * x + curve.a6()
    if legendre_symbol(rhs, p) == 1:
        return tonelli_shanks(rhs, p)
    else:
        return None

# Function to get small subgroups
def get_small_subgroups(E, threshold):
    Grp = E.abelian_group()
    g = Grp.gens()[0]
    numElements = g.order()
    factors = factor(numElements)
    small_subgroups = [factor[0] for factor in factors if factor[0] < threshold]
    return small_subgroups, numElements, factors

# Function to find and verify partial private keys
def find_and_verify_partial_private_keys(E, Q, small_subgroups):
    partial_private_keys = []
    for subgroup_order in small_subgroups:
        h = E.order() // subgroup_order
        print(f"Subgroup order: {subgroup_order}, h: {h}")
        try:
            hQ = h * E(Q)
            hg = h * E.gens()[0]
            print(f"hQ: {hQ}, hg: {hg}")
            partial_private_key = discrete_log(hQ, hg, operation='+')
            # Verify the partial private key on the twisted curve
            if hQ == partial_private_key * hg:
                partial_private_keys.append((partial_private_key, subgroup_order))
                print(f"Partial private key {partial_private_key} for subgroup order {subgroup_order} verified.")
            else:
                print(f"Partial private key {partial_private_key} for subgroup order {subgroup_order} failed verification.")
        except (ValueError, RuntimeError) as e:
            print(f"Error in discrete log computation: {e}")
            continue
    return partial_private_keys

# Function to check pairwise coprimeness
def check_pairwise_coprime(orders):
    for i in range(len(orders)):
        for j in range(i + 1, len(orders)):
            if gcd(orders[i], orders[j]) != 1:
                return False
    return True

# Define a threshold for small subgroups
threshold = 10**21

# List to store all verified partial private keys and their corresponding orders
all_partial_keys = []

# Iterate through the curves, get small subgroups, and find partial private keys
for name, curve in curves.items():
    small_subgroups, numElements, factors = get_small_subgroups(curve, threshold)
    print("-------------------------------------------------")
    print(f"{name}:")
    print(f"Order: {numElements}")
    print(f"Sub Groups: {factors}")
    print(f"Small subgroups: {small_subgroups}")
    print("-------------------------------------------------\n")

    # Use the provided x-coordinate and find a valid y-coordinate on the twisted curve
    Qx = pubx
    Qy = find_valid_y_on_twist(Qx, curve)
    if Qy is None:
        print(f"Could not find a valid y-coordinate on {name}")
        continue
    Q = (Qx, Qy)
    print(f"Using point Q = {Q} on {name}")

    # Find and verify partial private keys for the current curve
    partial_private_keys = find_and_verify_partial_private_keys(curve, Q, small_subgroups)
    print(f"Partial private keys for {name}")
    print(partial_private_keys)
    print("-------------------------------------------------\n")

    # Collect all partial private keys and orders for later combination
    all_partial_keys.extend(partial_private_keys)

# Filter out zero keys and ensure they are pairwise coprime
all_partial_keys = [(key, order) for key, order in all_partial_keys if key != 0]


