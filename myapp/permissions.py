from django.contrib.auth.decorators import user_passes_test
 
# Verifica si el usuario es administrador
def es_administrador(user):
    return user.is_authenticated and user.rol.nombre == 'Administrador'
 
# Decorador para vistas protegidas
def administrador_required(view_func):
    decorated_view_func = user_passes_test(es_administrador)(view_func)
    return decorated_view_func
 
from django.shortcuts import redirect
from django.contrib.auth.decorators import login_required
from functools import wraps
 
def role_required(role_name):
    def decorator(view_func):
        @wraps(view_func)
        @login_required
        def _wrapped_view(request, *args, **kwargs):
            if request.user.is_authenticated and request.user.rol and request.user.rol.nombre == role_name:
                return view_func(request, *args, **kwargs)
            else:
                # Aquí puedes redirigir a una página de acceso denegado o al home
                return redirect('accesoDenegado')  # O cualquier página que tú definas
        return _wrapped_view
    return decorator