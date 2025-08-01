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

    path('recurso/<int:recurso_id>/editar/video/', views.editar_video, name='editar_video'),
    path('recurso/<int:recurso_id>/editar/texto/', views.editar_texto, name='editar_texto'),


    path('usAdmin/detalleUsuarios', views.detalleUsuarios, name="detalleUsuarios-adm"),
    path('usuarios/<int:user_id>/editar/', views.editar_usuario, name='editar_usuario'),
    path('usuarios/<int:user_id>/eliminar/', views.eliminar_usuario, name='eliminar_usuario'),

    path('usAdmin/<int:id>/resetearContra', views.resetContra_usuario, name="resetear_contrasena"),
    path('recurso/<int:recurso_id>/', views.detalle_recurso, name='detalle_recurso'),
    path('unirse-curso/', views.unirse_curso, name='unirse_curso'),
    path('mis-cursos/', views.listar_cursos, name='listar_cursos'),
    path('detalleCursoEstudiante/<int:curso_id>/', views.detalle_cursoEstudiante, name='detalle_cursoEstudiante'),
    path('practica/<int:practica_id>/resolver/', views.resolver_practica, name='resolver_practica'),
    path('practica/<int:practica_id>/historial/', views.historial_intentos, name='historial_intentos'),


    path('curso/<int:curso_id>/modulo/nuevo/', views.crear_modulo, name='crear_modulo'),
    path('seccion/<int:seccion_id>/crear-recurso/', views.crear_recurso, name='crear_recurso'),
    path('practica/<int:recurso_id>/ver/', views.ver_practica, name='ver_practica'),
    path('practica/<int:recurso_id>/editar/', views.editar_practica, name='editar_practica'),
    path('biblioteca/practicas/', views.biblioteca_practicas, name='biblioteca_practicas'),

    path('competencia/<int:competencia_id>/panel/', views.panel_competencia, name='panel_competencia'),
    path('competencia/<int:competencia_id>/agregar-pregunta/', views.agregar_pregunta_competencia, name='agregar_pregunta_competencia'),

    path('competencia/<int:competencia_id>/iniciar/', views.iniciar_competencia, name='iniciar_competencia'),
    path('competencia/pregunta/<int:pregunta_id>/enviar/', views.enviar_pregunta, name='enviar_pregunta'),



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