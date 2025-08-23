from django.urls import path
from myapp import views

urlpatterns = [
    path('', views.inicio, name='inicio'),
    path('login/', views.vistaLogin, name='inicio'),
    path("logout/", views.signout, name="logout"),
    path('registro/', views.registroUser),
    path('acceso_denegado/', views.acceso_denegado, name="acceso_denegado"),
    # proyecto/urls.py
    path('editar-perfil/', views.editar_perfil, name='editar_perfil'),
    path('coming_soon/', views.comming_soon, name='comming_soon'),


    #path('inicio/',views.custom_login),
    path('inicio/', views.custom_login, name='custom_login'), #cambio de Paula
    path('usAdmin/', views.inicioAdmin, name="dashboard-adm"),
    path('usAdmin/detalleCursos', views.vistaAllCursos, name="detalleCursos-adm"),
    path('detalleCurso/<int:curso_id>/', views.detalle_curso, name='detalle_curso'),

    path('recurso/<int:recurso_id>/editar/video/', views.editar_video, name='editar_video'),
    path('recurso/<int:recurso_id>/editar/texto/', views.editar_texto, name='editar_texto'),


    path('usAdmin/crearUsuario/', views.crear_usuario_admin, name='crear_usuario_admin'),
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

    # Dispatcher / edición de práctica (decide a qué UI enviar)
    path('practicas/<int:recurso_id>/editar/', views.editar_practica, name='editar_practica'),

    # UIs específicas
    path('practicas/<int:recurso_id>/editar/jsme/', views.practica_jsme, name='practica_jsme'),
    path('practicas/<int:recurso_id>/editar/builder2d/', views.practica_builder2d, name='practica_builder2d'),


    # --- constructor 2D ---
    path('practicas/<int:recurso_id>/builder2d/', views.practica_builder2d, name='practica_builder2d'),
    path('practicas/<int:recurso_id>/guardar-builder2d/', views.guardar_practica_builder2d, name='guardar_practica_builder2d'),
    path('api/validar-armado/<int:recurso_id>/', views.validar_molecula_armada, name='validar_molecula_armada'),    
    path('api/resolve-smiles/', views.resolve_smiles, name='resolve_smiles'),
    


    # API de validación para el Constructor 2D
    #path('practicas/<int:recurso_id>/validar-armado/', views.validar_molecula_armada, name='validar_molecula_armada'),
    path('api/smiles-a-grafo/', views.smiles_a_grafo, name='smiles_a_grafo'),


    path('curso/<int:curso_id>/modulo/nuevo/', views.crear_modulo, name='crear_modulo'),
    path('seccion/<int:seccion_id>/crear-recurso/', views.crear_recurso, name='crear_recurso'),
    path('practica/<int:recurso_id>/ver/', views.ver_practica, name='ver_practica'),
    path('practica/<int:recurso_id>/editar/', views.editar_practica, name='editar_practica'),
    path('biblioteca/practicas/', views.biblioteca_practicas, name='biblioteca_practicas'),
    path('recurso/<int:recurso_id>/editar/html/', views.editar_html, name='editar_html'),


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
    path('registro-docente/', views.registroDocente, name='registro_docente'),#BY TOMMY
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
    path('estudiante/finalizar-cuestionario/<int:intento_id>/', views.finalizar_cuestionario_estudiante, name='finalizar_cuestionario_estudiante'),
    path('estudiante/resultado-cuestionario/<int:intento_id>/', views.resultado_cuestionario, name='resultado_cuestionario'),
    path('estudiante/historial-cuestionario/<int:cuestionario_id>/', views.historial_cuestionario_especifico, name='historial_cuestionario_especifico'),


    ######## ===== SISTEMA DE COMPETENCIAS - DOCENTE =====
    # Creación y edición básica
    path('docente/crear_competencia/', views.crear_competencia, name='crear_competencia'),
    path('docente/crear_competencia/<int:seccion_id>/', views.crear_competencia, name='crear_competencia_seccion'),
    path('docente/editar_competencia/<int:competencia_id>/', views.editar_competencia, name='editar_competencia'),
    path('docente/editar_competencia_recurso/<int:recurso_id>/', views.editar_competencia_por_recurso, name='editar_competencia_por_recurso'),
    
    # Gestión de datos de competencias
    path('docente/competencia/guardar-competencia/', views.guardar_competencia, name='guardar_competencia'),
    path('docente/competencia/actualizar-configuracion-competencia/', views.actualizar_configuracion_competencia, name='actualizar_configuracion_competencia'),
    
    # Gestión de preguntas
    path('docente/competencia/agregar-pregunta-competencia/', views.agregar_pregunta_competencia, name='agregar_pregunta_competencia'),
    path('docente/competencia/actualizar-pregunta-competencia/', views.actualizar_pregunta_competencia, name='actualizar_pregunta_competencia'),
    path('docente/competencia/eliminar-pregunta-competencia/<int:pregunta_id>/', views.eliminar_pregunta_competencia, name='eliminar_pregunta_competencia'),
    path('docente/competencia/duplicar-pregunta-competencia/', views.duplicar_pregunta_competencia, name='duplicar_pregunta_competencia'),
    path('docente/competencia/cambiar-tipo-pregunta-competencia/', views.cambiar_tipo_pregunta_competencia, name='cambiar_tipo_pregunta_competencia'),
    
    # Gestión de opciones de preguntas
    path('docente/competencia/agregar-opcion-competencia/', views.agregar_opcion_competencia, name='agregar_opcion_competencia'),
    path('docente/competencia/actualizar-opcion-competencia/', views.actualizar_opcion_competencia, name='actualizar_opcion_competencia'),
    path('docente/competencia/eliminar-opcion-competencia/<int:opcion_id>/', views.eliminar_opcion_competencia, name='eliminar_opcion_competencia'),
    
    # Control de competencia en vivo
    path('docente/competencia/abrir-sala-espera/<int:competencia_id>/', views.abrir_sala_espera, name='abrir_sala_espera'),
    path('docente/competencia/panel/<int:competencia_id>/', views.panel_control_competencia, name='panel_control_competencia'),
    path('docente/competencia/iniciar/<int:competencia_id>/', views.iniciar_competencia_live, name='iniciar_competencia_live'),
    path('docente/competencia/finalizar/<int:competencia_id>/', views.finalizar_competencia_live, name='finalizar_competencia_live'),
    path('docente/competencia/pausar/<int:competencia_id>/', views.pausar_competencia, name='pausar_competencia'),
    path('docente/competencia/reanudar/<int:competencia_id>/', views.reanudar_competencia, name='reanudar_competencia'),
    path('docente/competencia/siguiente-pregunta/<int:competencia_id>/', views.siguiente_pregunta_competencia, name='siguiente_pregunta_competencia'),
    path('docente/competencia/pregunta-anterior/<int:competencia_id>/', views.pregunta_anterior_competencia, name='pregunta_anterior_competencia'),
    path('docente/competencia/ir-pregunta/<int:competencia_id>/<int:pregunta_orden>/', views.ir_a_pregunta_competencia, name='ir_a_pregunta_competencia'),
    
    # Gestión de grupos
    path('docente/competencia/gestionar-grupos/<int:competencia_id>/', views.gestionar_grupos, name='gestionar_grupos'),
    path('docente/competencia/formar-grupos/<int:competencia_id>/', views.formar_grupos_automatico, name='formar_grupos_automatico'),
    path('docente/competencia/crear-grupo-manual/<int:competencia_id>/', views.crear_grupo_manual, name='crear_grupo_manual'),
    path('docente/competencia/editar-grupo/<int:grupo_id>/', views.editar_grupo, name='editar_grupo'),
    path('docente/competencia/eliminar-grupo/<int:grupo_id>/', views.eliminar_grupo, name='eliminar_grupo'),
    path('docente/competencia/mover-estudiante-grupo/', views.mover_estudiante_grupo, name='mover_estudiante_grupo'),
    path('docente/competencia/asignar-estudiante-grupo/', views.asignar_estudiante_grupo, name='asignar_estudiante_grupo'),
    path('docente/competencia/remover-estudiante-grupo/<int:participacion_id>/', views.remover_estudiante_grupo, name='remover_estudiante_grupo'),
    
    # Monitoreo y resultados
    path('docente/competencia/participantes/<int:competencia_id>/', views.obtener_participantes_competencia, name='obtener_participantes_competencia'),
    path('docente/competencia/resultados/<int:competencia_id>/', views.resultados_competencia, name='resultados_competencia'),
    path('docente/competencia/resultados-detallados/<int:competencia_id>/', views.resultados_detallados_competencia, name='resultados_detallados_competencia'),
    path('docente/competencia/exportar/<int:competencia_id>/', views.exportar_resultados_competencia, name='exportar_resultados_competencia'),
    path('docente/competencia/exportar-excel/<int:competencia_id>/', views.exportar_resultados_excel, name='exportar_resultados_excel'),
    path('docente/competencia/estadisticas/<int:competencia_id>/', views.estadisticas_competencia, name='estadisticas_competencia'),
    
    # Gestión avanzada
    path('docente/competencia/duplicar/<int:competencia_id>/', views.duplicar_competencia, name='duplicar_competencia'),
    path('docente/competencia/reiniciar/<int:competencia_id>/', views.reiniciar_competencia, name='reiniciar_competencia'),
    path('docente/competencia/cancelar/<int:competencia_id>/', views.cancelar_competencia, name='cancelar_competencia'),
    path('docente/competencia/archivar/<int:competencia_id>/', views.archivar_competencia, name='archivar_competencia'),
    path('docente/competencia/restaurar/<int:competencia_id>/', views.restaurar_competencia, name='restaurar_competencia'),
    
    # Biblioteca y plantillas
    path('docente/biblioteca-competencias/', views.biblioteca_competencias, name='biblioteca_competencias'),
    path('docente/competencia/crear-desde-plantilla/<int:plantilla_id>/', views.crear_desde_plantilla, name='crear_competencia_desde_plantilla'),
    path('docente/competencia/guardar-como-plantilla/<int:competencia_id>/', views.guardar_como_plantilla, name='guardar_como_plantilla'),
    

    # ===== COMPETENCIAS - ESTUDIANTE =====
    
    # Acceso inicial
    path('competencia/acceder/<str:pin>/', views.acceder_competencia, name='acceder_competencia'),
    path('competencia/buscar-por-codigo/', views.buscar_competencia_codigo, name='buscar_competencia_codigo'),
    path('competencia/unirse/<int:competencia_id>/', views.unirse_competencia, name='unirse_competencia'),
    path('competencia/salir/<int:participacion_id>/', views.salir_competencia, name='salir_competencia'),
    
    # Gestión de grupos (estudiante)
    path('competencia/seleccionar-grupo/<int:competencia_id>/', views.seleccionar_grupo, name='seleccionar_grupo'),
    path('competencia/crear-grupo/<int:competencia_id>/', views.crear_grupo_estudiante, name='crear_grupo_estudiante'),
    path('competencia/unirse-grupo/<int:grupo_id>/', views.unirse_grupo, name='unirse_grupo'),
    path('competencia/salir-grupo/<int:participacion_id>/', views.salir_grupo, name='salir_grupo'),
    path('competencia/info-grupo/<int:grupo_id>/', views.info_grupo, name='info_grupo'),
    
    # Participación en competencia
    path('competencia/sala-espera/<int:participacion_id>/', views.sala_espera_competencia, name='sala_espera_competencia'),
    path('competencia/participar/<int:participacion_id>/', views.participar_competencia, name='participar_competencia'),
    path('competencia/responder/', views.responder_competencia, name='responder_competencia'),
    path('competencia/finalizar-participacion/<int:participacion_id>/', views.finalizar_participacion, name='finalizar_participacion'),
    path('competencia/abandonar/<int:participacion_id>/', views.abandonar_competencia, name='abandonar_competencia'),
    
    # Resultados y historial (estudiante)
    path('competencia/resultados-participacion/<int:participacion_id>/', views.resultados_participacion, name='resultados_participacion'),
    path('competencia/resultados-grupo/<int:grupo_id>/', views.resultados_grupo, name='resultados_grupo'),
    path('competencia/historial-competencias/', views.historial_competencias_estudiante, name='historial_competencias_estudiante'),
    path('competencia/certificado/<int:participacion_id>/', views.generar_certificado, name='generar_certificado'),
    
    # ===== AJAX Y TIEMPO REAL =====
    
    # Estado y sincronización
    path('competencia/estado/<int:competencia_id>/', views.obtener_estado_competencia, name='obtener_estado_competencia'),
    path('competencia/info-participacion/<int:participacion_id>/', views.obtener_info_participacion, name='obtener_info_participacion'),
    path('competencia/ping/<int:participacion_id>/', views.ping_participacion, name='ping_participacion'),
    path('competencia/heartbeat/<int:competencia_id>/', views.heartbeat_competencia, name='heartbeat_competencia'),
    
    # Rankings y estadísticas en tiempo real
    path('competencia/ranking/<int:competencia_id>/', views.obtener_ranking_tiempo_real, name='obtener_ranking_tiempo_real'),
    path('competencia/ranking-grupo/<int:competencia_id>/', views.obtener_ranking_grupos, name='obtener_ranking_grupos'),
    path('competencia/progreso/<int:competencia_id>/', views.obtener_progreso_competencia, name='obtener_progreso_competencia'),
    path('competencia/estadisticas-live/<int:competencia_id>/', views.estadisticas_live, name='estadisticas_live'),
    
    # Preguntas en tiempo real
    path('competencia/pregunta-actual/<int:competencia_id>/', views.obtener_pregunta_actual, name='obtener_pregunta_actual'),
    path('competencia/siguiente-pregunta-info/<int:competencia_id>/', views.info_siguiente_pregunta, name='info_siguiente_pregunta'),
    path('competencia/tiempo-restante/<int:competencia_id>/', views.tiempo_restante_competencia, name='tiempo_restante_competencia'),
    
    # Chat y comunicación (opcional)
    path('competencia/chat/<int:competencia_id>/', views.chat_competencia, name='chat_competencia'),
    path('competencia/enviar-mensaje/', views.enviar_mensaje_chat, name='enviar_mensaje_chat'),
    path('competencia/obtener-mensajes/<int:competencia_id>/', views.obtener_mensajes_chat, name='obtener_mensajes_chat'),
    
    # ===== APIs ADICIONALES =====
    
    # API para apps móviles
    path('api/competencia/listar/', views.api_listar_competencias, name='api_listar_competencias'),
    path('api/competencia/<int:competencia_id>/info/', views.api_info_competencia, name='api_info_competencia'),
    path('api/competencia/unirse/', views.api_unirse_competencia, name='api_unirse_competencia'),
    path('api/competencia/responder/', views.api_responder_competencia, name='api_responder_competencia'),
    
    # Búsqueda y filtros
    path('competencia/buscar/', views.buscar_competencias, name='buscar_competencias'),
    path('competencia/filtrar/', views.filtrar_competencias, name='filtrar_competencias'),
    
    # ===== WEBHOOKS Y NOTIFICACIONES =====
    
    # Notificaciones
    path('competencia/notificar-inicio/<int:competencia_id>/', views.notificar_inicio_competencia, name='notificar_inicio_competencia'),
    path('competencia/enviar-recordatorio/<int:competencia_id>/', views.enviar_recordatorio, name='enviar_recordatorio'),
    
    # Webhooks para integraciones
    path('webhook/competencia/inicio/', views.webhook_inicio_competencia, name='webhook_inicio_competencia'),
    path('webhook/competencia/finalizacion/', views.webhook_finalizacion_competencia, name='webhook_finalizacion_competencia'),
    


]