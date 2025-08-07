from django import template
import json

register = template.Library()

@register.filter
def split(value, arg):
    """Divide una cadena por el separador especificado"""
    if value:
        return str(value).split(str(arg))
    return []

@register.filter
def trim(value):
    """Elimina espacios en blanco de una cadena"""
    if value:
        return str(value).strip()
    return value

@register.filter
def get_item(dictionary, key):
    """Obtiene un elemento de un diccionario"""
    return dictionary.get(key)

@register.filter
def to_json(value):
    """Convierte un valor a JSON"""
    return json.dumps(value)

@register.filter
def make_list(value):
    """Convierte un rango en una lista"""
    try:
        return list(range(int(value)))
    except:
        return []

@register.filter
def add_num(value, arg):
    """Suma un número a otro"""
    try:
        return int(value) + int(arg)
    except:
        return value

@register.filter
def stringformat_char(value):
    """Convierte un número a su caracter ASCII correspondiente"""
    try:
        return chr(int(value))
    except:
        return value

@register.filter
def mul(value, arg):
    """Multiplica dos valores"""
    try:
        return float(value) * float(arg)
    except (TypeError, ValueError):
        return 0

@register.filter
def div(value, arg):
    """Divide dos valores"""
    try:
        if float(arg) == 0:
            return 0
        return float(value) / float(arg)
    except (TypeError, ValueError, ZeroDivisionError):
        return 0

@register.filter
def floatformat_zero(value):
    """Formatea un float sin decimales si es entero"""
    try:
        val = float(value)
        if val == int(val):
            return int(val)
        return round(val, 1)
    except (TypeError, ValueError):
        return value