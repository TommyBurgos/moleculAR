from django.urls import path
from myapp import views

urlpatterns = [
    path('', views.inicio, name='inicio'),
    path('login/', views.vistaLogin, name='inicio'),
    path("logout/", views.signout, name="logout"),
    path('registro/', views.registroUser),
    path('acceso_denegado/', views.acceso_denegado, name="acceso_denegado"),
    #path('inicio/',views.custom_login),
    path('inicio/', views.custom_login, name='custom_login'), #cambio de Paula
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


    #paula
    #Profesor
    path('docente/crear_cuestionario/', views.instructorQuiz, name='crearCuestionario'),
    path('docente/crear_cuestionario/<int:seccion_id>/', views.instructorQuiz, name='crearCuestionarioSeccion'),
    path('docente/editar_cuestionario/<int:cuestionario_id>/', views.editarCuestionario, name='editarCuestionario'),
    path('docente/editarCuestionario/<int:recurso_id>/', views.editar_cuestionario, name='editar_cuestionario'), #TOMMY> MI VERSION PARA QUE FUNCIONE AL AGREGAR RECURSO

    
    
    path('docente/cuestionario/guardar-cuestionario/', views.guardar_cuestionario, name='guardarCuestionario'),#NOTA TOMMY: ME SALE QUE NO FUNCIONA ETA RUTA.
    path('docente/cuestionario/actualizar-cuestionario/<int:cuestionario_id>/', views.actualizarCuestionarioExistente, name='actualizarCuestionarioExistente'),
    path('docente/cuestionario/actualizar-configuracion/', views.actualizarConfiguracion, name='actualizarConfiguracion'),#NOTA TOMMY: ME SALE QUE NO FUNCIONA ETA RUTA. al parecer por el required post.
    path('docente/cuestionario/finalizar/', views.finalizar_cuestionario, name='finalizarCuestionario'),
    
    path('docente/cuestionario/agregar-pregunta/', views.agregarPregunta, name='agregarPregunta'),
    path('docente/cuestionario/actualizar-pregunta/', views.actualizar_pregunta, name='actualizarPregunta'),
    path('docente/cuestionario/eliminar-pregunta/<int:pregunta_id>/', views.eliminarPregunta, name='eliminarPregunta'),
    path('docente/cuestionario/cambiar-tipo-pregunta/', views.cambiarTipoPregunta, name='cambiarTipoPregunta'),
    path('docente/cuestionario/duplicar-pregunta/', views.duplicar_pregunta, name='duplicarPregunta'),
    
    path('docente/cuestionario/agregar-opcion/', views.agregarOpcion, name='agregarOpcion'),
    path('docente/cuestionario/actualizar-opcion/', views.actualizarOpcion, name='actualizarOpcion'),
    path('docente/cuestionario/eliminar-opcion/<int:opcion_id>/', views.eliminarOpcion, name='eliminarOpcion'),
    
    path('docente/cuestionario/subir-recurso/', views.subir_recurso, name='subirRecurso'),
    path('docente/cuestionario/eliminar-recurso-pregunta/', views.eliminar_recurso_pregunta, name='eliminarRecursoPregunta'),
    path('docente/cuestionario/crear-tipos/', views.crearTiposPregunta, name='crearTiposPregunta'),
    path('docente/cuestionario/inicializar-tipos/', views.inicializar_tipos_pregunta, name='inicializarTiposPregunta'),
    path('docente/biblioteca-cuestionarios/', views.biblioteca_cuestionarios, name='bibliotecaCuestionarios'),
    path('docente/asignar-cuestionario/<int:cuestionario_id>/', views.asignar_cuestionario_seccion, name='asignarCuestionarioSeccion'),

    #Estudiante
    path('estudiante/biblioteca-cuestionarios/', views.biblioteca_cuestionarios_estudiante, name='biblioteca_cuestionarios_estudiante'),
    path('estudiante/iniciar-cuestionario/<int:cuestionario_id>/', views.iniciar_cuestionario, name='iniciar_cuestionario'),
    path('estudiante/realizar-cuestionario/<int:intento_id>/', views.realizar_cuestionario, name='realizar_cuestionario'),
    path('estudiante/guardar-respuesta/', views.guardar_respuesta, name='guardar_respuesta'),
    path('estudiante/finalizar-cuestionario/<int:intento_id>/', views.finalizar_cuestionario, name='finalizar_cuestionario'),
    path('estudiante/resultado-cuestionario/<int:intento_id>/', views.resultado_cuestionario, name='resultado_cuestionario'),
    path('estudiante/historial-cuestionarios/', views.historial_cuestionarios, name='historial_cuestionarios'),


]