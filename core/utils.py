from math import floor, log10
from numpy import isnan
from .to_precision import std_notation, sci_notation

# Contains functions for rounding values

superscript = {
        '0': '⁰',
        '1': '¹',
        '2': '²',
        '3': '³',
        '4': '⁴',
        '5': '⁵',
        '6': '⁶',
        '7': '⁷',
        '8': '⁸',
        '9': '⁹',
        ',': '˒',
        'e': 'ᵉ',
        'r': 'ʳ',
        '-': '⁻',
        '.': 'ˑ'
}


def superscriptize_string(cdc):
    atomisé = tuple(superscript[c] for c in cdc)
    return "".join(atomisé)


def superscriptize_number(nombre, decimal_comma):
    cdc = str(nombre)
    cdc = cdc.replace(".", ",") if decimal_comma else cdc
    return superscriptize_string(cdc)

def round_with_sense(v, digit, decimal_comma=False):
    # retourne une chaîne de caractère.
    if isnan(v):
        return ("/")
    if digit is None:
        strv = str(v)
        return strv.replace(".", ",") if decimal_comma else strv
    s = std_notation(v, digit)
    if len(s) > digit+3:
        s = sci_notation(v, digit, delimiter='e')
        b, p = s.split("e")
        b = b.replace(".", ",") if decimal_comma else b
        return b+"·10"+superscriptize_number(p, decimal_comma)
    else:
        if s[-1] == ".":
            s = s[:-1]
        return s.replace(".", ",") if decimal_comma else s

def round_according_to_std(value, std, decimal_comma):    
    # Arrondi `value` en fonction de l'incertitude selon le GUM
    if isnan(value):
        return ("/")
    if isnan(std):
        digits = 2
    elif std == 0:
        digits = 5
    else:
        digits = min(max(n_digits_according_to_std(value, std), 2), 5)
    return round_with_sense(value, digits, decimal_comma)


def n_digits_according_to_std(val, std):
    # Retourn le nb de chiffres significatifs pour val
    # en fonction de l'incertitude selon le GUM
    return int(floor(log10(abs(val)))) - int(floor(log10(abs(2*std)))) + 2


