from django.shortcuts import redirect
from django.contrib.auth.decorators import user_passes_test, login_required
from functools import wraps

# Verifica si el usuario es administrador
def es_administrador(user):
    print("Dentro de la funcion es_administrador: ")
    print(user.is_authenticated and user.rol.nombre == 'Administrador')
    return user.is_authenticated and user.rol.nombre == 'Administrador'

# Decorador para vistas protegidas (solo administrador)
def administrador_required(view_func):
    decorated_view_func = user_passes_test(es_administrador)(view_func)
    return decorated_view_func

# Decorador para uno o más roles permitidos
def role_required(role_names):
    if isinstance(role_names, str):
        print(role_names)        
        role_names = [role_names]
        print(f"el rol name ahora es: {role_names}")

    def decorator(view_func):
        print("Entré al decorator")
        @wraps(view_func)        
        def _wrapped_view(request, *args, **kwargs):
            print(f"el rol name ahora es: {role_names}")
            print(request.user.is_authenticated and request.user.rol.nombre in role_names)
            print(request.user.is_authenticated)            
            if request.user.is_authenticated and request.user.rol.nombre in role_names:
                return view_func(request, *args, **kwargs)
            return redirect('acceso_denegado') 
        return _wrapped_view
    return decorator
