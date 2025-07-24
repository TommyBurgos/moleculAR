from django.urls import path
from myapp import views

urlpatterns = [
    path('', views.inicio, name='inicio'),
    path('login/', views.vistaLogin, name='inicio'),
    path("logout/", views.signout, name="logout"),
    path('registro/', views.registroUser),
    path('inicio/',views.custom_login),
    path('usAdmin/', views.inicioAdmin, name="dashboard-adm"),
    path('usAdmin/detalleCursos', views.vistaAllCursos, name="detalleCursos-adm"),
    path('detalleCurso/<int:curso_id>/', views.detalle_curso, name='detalle_curso'),
    path('recurso/<int:recurso_id>/', views.detalle_recurso, name='detalle_recurso'),

    path('curso/<int:curso_id>/modulo/nuevo/', views.crear_modulo, name='crear_modulo'),
    path('seccion/<int:seccion_id>/crear-recurso/', views.crear_recurso, name='crear_recurso'),
    path('practica/<int:recurso_id>/ver/', views.ver_practica, name='ver_practica'),
    path('practica/<int:recurso_id>/editar/', views.editar_practica, name='editar_practica'),




    path('usAdmin/crearCurso', views.vistaCrearCurso, name="crearCurso-adm"),
    path("obtener-presigned-url/", views.obtener_presigned_url, name="obtener_presigned_url"),

    path('usEstudiante/', views.inicioEstudiante, name="student_dashboard"),
    path('usDocente/', views.inicioDocente, name="teacher_dashboard"),
    path("modificar-molecula/", views.modificar_molecula, name="modificar_molecula"),
    path("editor/", views.vista_modificador, name="vista_modificador"),
    path('editor-dibujo/', views.editor_dibujo, name='editor_dibujo'),
    path('procesar-molecula/', views.procesar_molecula, name='procesar_molecula'),
    path('registroExitoso/',views.registroExitoso),


]