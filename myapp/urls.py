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
    path('curso/<int:curso_id>/modulo/nuevo/', views.crear_modulo, name='crear_modulo'),
    path('seccion/<int:seccion_id>/crear-recurso/', views.crear_recurso, name='crear_recurso'),


    path('usAdmin/crearCurso', views.vistaCrearCurso, name="crearCurso-adm"),
    path("obtener-presigned-url/", views.obtener_presigned_url, name="obtener_presigned_url"),

    path('usEstudiante/', views.inicioEstudiante, name="student_dashboard"),
    path('usDocente/', views.inicioDocente, name="teacher_dashboard"),
    path("modificar-molecula/", views.modificar_molecula, name="modificar_molecula"),
    path("editor/", views.vista_modificador, name="vista_modificador"),
    path('editor-dibujo/', views.editor_dibujo, name='editor_dibujo'),
    path('procesar-molecula/', views.procesar_molecula, name='procesar_molecula'),
    path('registroExitoso/',views.registroExitoso),

    #paula
    path('docente/crear_cuestionario/', views.instructorQuiz, name='crearCuestionario'),
    
    # URLs para manejar cuestionarios con sección específica
    path('docente/crear_cuestionario/<int:seccion_id>/', views.instructorQuiz, name='crearCuestionarioSeccion'),
    
    # Editar cuestionario existente
    path('docente/editar_cuestionario/<int:cuestionario_id>/', views.editarCuestionario, name='editarCuestionario'),
    
    # APIs AJAX para las funcionalidades del cuestionario
    path('docente/cuestionario/agregar-pregunta/', views.agregarPregunta, name='agregarPregunta'),
    path('docente/cuestionario/agregar-opcion/', views.agregarOpcion, name='agregarOpcion'),
    path('docente/cuestionario/actualizar-opcion/', views.actualizarOpcion, name='actualizarOpcion'),
    path('docente/cuestionario/actualizar-pregunta/', views.actualizarPregunta, name='actualizarPregunta'),
    path('docente/cuestionario/actualizar-configuracion/', views.actualizarConfiguracion, name='actualizarConfiguracion'),
    path('docente/cuestionario/eliminar-opcion/<int:opcion_id>/', views.eliminarOpcion, name='eliminarOpcion'),
    path('docente/cuestionario/eliminar-pregunta/<int:pregunta_id>/', views.eliminarPregunta, name='eliminarPregunta'),
    path('docente/cuestionario/crear-tipos/', views.crearTiposPregunta, name='crearTiposPregunta'),
    # URL para cambiar tipo de pregunta
    path('docente/cuestionario/cambiar-tipo-pregunta/', views.cambiarTipoPregunta, name='cambiarTipoPregunta'),
    # URL para guardar cuestionario completo
    path('docente/cuestionario/guardar-cuestionario/', views.guardarCuestionarioCompleto, name='guardarCuestionarioCompleto'),


]