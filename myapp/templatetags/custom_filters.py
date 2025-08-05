from django import template
import json

register = template.Library()

@register.filter
def split(value, arg):
    """Divide una cadena por el separador especificado"""
    if value:
        return value.split(arg)
    return []

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