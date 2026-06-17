from math import sin, cos

def get_portal_params(mx, mh2, gX, beta, mh1=125.0, vEW=246.0):

    vPhi = 2.0 * mx / gX
    
    c = cos(beta)
    s = sin(beta)

    lambdaH   = (mh1**2 * c**2 + mh2**2 * s**2) / (2.0 * vEW**2)
    lambdaPhi = (mh1**2 * s**2 + mh2**2 * c**2) / (2.0 * vPhi**2)
    lambdaM   = (mh2**2 - mh1**2) * s * c / (vEW * vPhi)
    
    mu    = - lambdaH * vEW**2 - 0.5 * lambdaM * vPhi**2
    muPhi = - lambdaPhi * vPhi**2 - 0.5 * lambdaM * vEW**2

    return {
        "vPhi": float(vPhi),
        "lambdaH": float(lambdaH),
        "lambdaPhi": float(lambdaPhi),
        "lambdaM": float(lambdaM),
        "mu": float(mu),
        "muPhi": float(muPhi),
        "mh1": float(mh1),
        "mh2": float(mh2),
        "vEW": float(vEW),
        "beta": float(beta),
        "gX": float(gX),
        "mx": float(mx),
    }
