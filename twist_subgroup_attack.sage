from sage.all import *

# Define the prime modulus
p = 115792089237316195423570985008687907853269984665640564039457584007908834671663
Fp = GF(p)  # Define the field over the prime modulus

# Define the elliptic curves
curves = {
    "Original Curve": EllipticCurve(Fp, [0, 7]),
    "Twist Curve 1": EllipticCurve(Fp, [0, -9]),
    "Twist Curve 2": EllipticCurve(Fp, [0, -7]),
    "Twist Curve 3": EllipticCurve(Fp, [0, -6]),
    "Twist Curve 4": EllipticCurve(Fp, [0, -1]),
    "Twist Curve 5": EllipticCurve(Fp, [0, 1]),
    "Twist Curve 6": EllipticCurve(Fp, [0, 4]),
    "Twist Curve 7": EllipticCurve(Fp, [0, 5]),
    "Twist Curve 8": EllipticCurve(Fp, [0, 6]),
    "Twist Curve 9": EllipticCurve(Fp, [0, 8]),
    "Twist Curve 10": EllipticCurve(Fp, [0, 9]),
    "Twist Curve 11": EllipticCurve(Fp, [0, 12]),
    "Twist Curve 12": EllipticCurve(Fp, [0, 14]),
    "Twist Curve 13": EllipticCurve(Fp, [0, 17])
}

# Convert public key x-coordinate from hex. Input your public key X Coordinate
pubx_int = int("1abc12334509fc", 16)

def validate_curve_points():
    x = Fp(pubx_int)
    for name, curve in curves.items():
        print(f"\nValidating points on {name}:")
        print(f"  Discriminant: {curve.discriminant()}")
        rhs = x^3 + curve.a6()
        try:
            ys = rhs.sqrt(all=True)
            for y in ys:
                point = (x, y)
                if curve.is_on_curve(x, y):
                    print(f"  Valid point on {name}: {point}")
                else:
                    print(f"  Invalid point on {name}: {point}")
        except (ValueError, ArithmeticError) as e:
            print(f"  No valid points found for x = {x} on {name}. Error: {e}")

def process_curve_groups():
    for name, curve in curves.items():
        print(f"\nProcessing group details for {name}:")
        Grp = curve.abelian_group()
        g = Grp.gens()[0]
        numElements = g.order()
        factorization = factor(numElements)
        print(f"  Group order: {numElements}")
        print(f"  Factors: {factorization}")
        for fac, exp in factorization:
            n = fac**exp
            h = numElements // n
            subgroup_gen = h * g
            print(f"  Subgroup order: {n}, Generator: {subgroup_gen}")
            print(f"  Order of subgroup generator: {subgroup_gen.order()}")

# Validate points on all twist curves
validate_curve_points()

# Process all curve groups
process_curve_groups()
