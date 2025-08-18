from django.db import models
from django.contrib.auth.models import AbstractUser
from django.contrib.contenttypes.models import ContentType
from django.contrib.contenttypes.fields import GenericForeignKey
import string
import random
from django.utils import timezone

# Create your models here.
class Rol(models.Model):
    nombre = models.CharField(max_length=50, unique=True)
    descripcion = models.TextField(blank=True)

    def __str__(self):
        return self.nombre

class User(AbstractUser):
    imgPerfil = models.ImageField(default='imageDefault.png', upload_to='users/')    
    rol = models.ForeignKey('Rol', on_delete=models.CASCADE, null=True, blank=True)

class CategoriaCurso(models.Model):
    nombre = models.CharField(max_length=100, unique=True)
    descripcion = models.TextField(blank=True)

    def __str__(self):
        return self.nombre

import string, random

def generar_codigo_aleatorio():
    return ''.join(random.choices(string.ascii_uppercase + string.digits, k=10))

class Curso(models.Model):
    titulo = models.CharField(max_length=200)
    descripcion = models.TextField()
    imagen_portada = models.ImageField(upload_to='cursos/', blank=True, null=True)
    profesor = models.ForeignKey('User', on_delete=models.CASCADE, related_name='cursos_creados')
    categoria = models.ForeignKey('CategoriaCurso', on_delete=models.SET_NULL, null=True, blank=True)
    fecha_creacion = models.DateTimeField(auto_now_add=True)
    codigo_acceso = models.CharField(max_length=10, unique=True, blank=True, null=True)

    def save(self, *args, **kwargs):
        if not self.codigo_acceso:
            self.codigo_acceso = generar_codigo_aleatorio()
        super().save(*args, **kwargs)

    def __str__(self):
        return self.titulo

class InscripcionCurso(models.Model):
    estudiante = models.ForeignKey('User', on_delete=models.CASCADE)
    curso = models.ForeignKey('Curso', on_delete=models.CASCADE)
    fecha_inscripcion = models.DateTimeField(auto_now_add=True)

    class Meta:
        unique_together = ('estudiante', 'curso')  # Para evitar inscripciones duplicadas


class Seccion(models.Model):
    curso = models.ForeignKey('Curso', on_delete=models.CASCADE, related_name='secciones')
    titulo = models.CharField(max_length=150)
    descripcion = models.TextField(blank=True)
    orden = models.PositiveIntegerField(default=1)  # Para ordenar secciones

    def __str__(self):
        return f"{self.titulo} (Curso: {self.curso.titulo})"


class TipoRecurso(models.Model):
    nombre = models.CharField(max_length=50, unique=True)  # Ej: video, texto, cuestionario, simulacion
    descripcion = models.TextField(blank=True)

    def __str__(self):
        return self.nombre


class Recurso(models.Model):
    #seccion = models.ForeignKey('Seccion', on_delete=models.CASCADE, related_name='recursos')
    seccion = models.ForeignKey('Seccion', on_delete=models.CASCADE, related_name='recursos', 
                               null=True, blank=True) #ESTA SOLO ES DE PRUEBA, EN REALIDAD SE USA LA LINEA DE ARRIBA
    tipo = models.ForeignKey('TipoRecurso', on_delete=models.PROTECT)
    titulo = models.CharField(max_length=150)
    descripcion = models.TextField(blank=True)
    archivo = models.FileField(upload_to='recursos/', blank=True, null=True)
    imagen = models.ImageField(upload_to='recursos/imgs/', blank=True, null=True)
    orden = models.PositiveIntegerField(default=1)  # Para ordenar recursos dentro de la sección
    fecha_creacion = models.DateTimeField(auto_now_add=True)
    contenido_texto = models.TextField(blank=True, null=True)  #NUEVO CAMPO para recurso tipo texto
    visible_biblioteca = models.BooleanField(default=False)  # Marco si va a la biblioteca pública
    video_url = models.URLField(blank=True, null=True)
    contenido_html = models.TextField(blank=True, null=True)  # Código HTML puro


    creado_por = models.ForeignKey('User', on_delete=models.SET_NULL, null=True, blank=True, 
                                   related_name='recursos_creados')


    def __str__(self):
        return f"{self.titulo} ({self.tipo.nombre})"

class Cuestionario(models.Model):
    recurso = models.OneToOneField('Recurso', on_delete=models.CASCADE, related_name='cuestionario')
    instrucciones = models.TextField(blank=True)
    tiempo_limite = models.PositiveIntegerField(help_text='Tiempo en minutos', default=10)
    
    # NUEVOS CAMPOS A AGREGAR:
    fecha_apertura = models.DateTimeField(null=True, blank=True, help_text='Fecha y hora de apertura del cuestionario')
    fecha_cierre = models.DateTimeField(null=True, blank=True, help_text='Fecha y hora de cierre del cuestionario')
    puntaje_total = models.DecimalField(max_digits=6, decimal_places=2, default=0.00, help_text='Puntaje total del cuestionario')
    calificacion_automatica = models.BooleanField(default=True, help_text='Si se califica automáticamente o requiere revisión manual')
    intentos_permitidos = models.PositiveIntegerField(default=1, help_text='Número de intentos permitidos')
    mostrar_resultados = models.BooleanField(default=True, help_text='Mostrar resultados al estudiante después de completar')
    orden_aleatorio = models.BooleanField(default=False, help_text='Mostrar preguntas en orden aleatorio')

    def __str__(self):
        return f"Cuestionario - {self.recurso.titulo}"
    
    def calcular_puntaje_total(self):
        "Calcula el puntaje total sumando todas las preguntas"
        return self.preguntas.aggregate(
            total=models.Sum('puntaje')
        )['total'] or 0
    
    def tiene_preguntas_manuales(self):
        "Verifica si tiene preguntas que requieren calificación manual"
        return self.preguntas.filter(
            tipo__nombre__in=['respuesta_abierta', 'simulador_2d', 'simulador_3d']
        ).exists()
    
    def esta_disponible(self):
        "Verifica si el cuestionario está disponible según las fechas"
        from django.utils import timezone
        ahora = timezone.now()
        
        if self.fecha_apertura and ahora < self.fecha_apertura:
            return False
        if self.fecha_cierre and ahora > self.fecha_cierre:
            return False
        return True

class Practica(models.Model):
    recurso = models.OneToOneField('Recurso', on_delete=models.CASCADE, related_name='practica')
    titulo_objetivo = models.CharField(max_length=150)
    modelo_objetivo = models.ForeignKey('Modelo', on_delete=models.SET_NULL, null=True, blank=True)
    instrucciones = models.TextField(blank=True)

    def __str__(self):
        return f"Práctica - {self.recurso.titulo}"

import json
from typing import Union, Optional

class PracticaConfig(models.Model):
    TIPO_PRACTICA = [('jsme','Editor JSME'), ('builder2d','Constructor 2D')]
    recurso = models.OneToOneField(Recurso, on_delete=models.CASCADE, related_name='practica_cfg')
    practica_tipo = models.CharField(max_length=20, choices=TIPO_PRACTICA, default='jsme')
    objetivo_smiles = models.TextField(blank=True, null=True)
    modelo = models.ForeignKey('Modelo', blank=True, null=True, on_delete=models.SET_NULL)
    # antes: JSONField
    inventario_json = models.TextField(blank=True, null=True)

    # helpers (opcionales)
    def get_inventario(self):
        try:
            return json.loads(self.inventario_json) if self.inventario_json else None
        except Exception:
            return None    
    def set_inventario(self, data: Optional[Union[dict, list]]):

        self.inventario_json = json.dumps(data) if data is not None else None

class IntentoArmado(models.Model):
    practica = models.ForeignKey('Recurso', on_delete=models.CASCADE, related_name='intentos_armado')
    estudiante = models.ForeignKey('User', on_delete=models.SET_NULL, null=True, blank=True)
    # antes: JSONField
    graph_json = models.TextField()   # guardaremos el JSON como string
    smiles_generado = models.TextField(blank=True)
    es_correcto = models.BooleanField(default=False)
    creado_en = models.DateTimeField(auto_now_add=True)

class Modelo(models.Model):
    titulo = models.CharField(max_length=150)
    descripcion = models.TextField(blank=True)
    tipo = models.CharField(max_length=10, choices=[('2D', '2D'), ('3D', '3D')])
    archivo = models.FileField(upload_to='modelos/')
    smiles = models.TextField(blank=True)  # Si aplica
    creado_por = models.ForeignKey('User', on_delete=models.SET_NULL, null=True)
    visible_biblioteca = models.BooleanField(default=False)
    fecha_creacion = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return self.titulo


class MoleculaEstudiante(models.Model):
    estudiante = models.ForeignKey('User', on_delete=models.CASCADE)
    practica = models.ForeignKey('Practica', on_delete=models.CASCADE)
    smiles = models.TextField()
    archivo = models.FileField(upload_to='moleculas/', blank=True)  # si quieres guardar archivo .mol o .cml
    fecha_envio = models.DateTimeField(auto_now_add=True)
    es_correcta = models.BooleanField(null=True, blank=True)

    def __str__(self):
        return f"Molecula de {self.estudiante.username} para {self.practica.recurso.titulo}"


#Competencias
def generar_pin_competencia():
    """Genera un PIN único de 6 dígitos para la competencia"""
    return ''.join(random.choices(string.digits, k=6))

class Competencia(models.Model):
    # Relación con recurso (igual que antes)
    recurso = models.OneToOneField('Recurso', on_delete=models.CASCADE, related_name='competencia')
    
    # Configuración básica
    instrucciones = models.TextField(blank=True, help_text='Instrucciones para los participantes')
    codigo_acceso = models.CharField(max_length=10, unique=True, blank=True, help_text='Código para unirse a la competencia')
    pin_acceso = models.CharField(max_length=6, unique=True, blank=True, help_text='PIN de acceso para estudiantes')
    
    # Configuración de modalidad
    MODALIDAD_CHOICES = [
        ('individual', 'Individual'),
        ('grupal', 'Grupal'),
    ]
    modalidad = models.CharField(max_length=10, choices=MODALIDAD_CHOICES, default='individual')
    
    # Configuración para modalidad grupal
    max_miembros_grupo = models.PositiveIntegerField(default=4, help_text='Máximo de miembros por grupo')
    grupos_aleatorios = models.BooleanField(default=False, help_text='Formar grupos automáticamente de forma aleatoria')
    grupos_abiertos = models.BooleanField(default=True, help_text='Permitir que estudiantes formen sus propios grupos')
    
    # Estados de la competencia
    ESTADO_CHOICES = [
        ('configuracion', 'En Configuración'),
        ('esperando', 'Esperando Participantes'),
        ('activa', 'Activa'),
        ('finalizada', 'Finalizada'),
        ('cancelada', 'Cancelada'),
    ]
    estado = models.CharField(max_length=15, choices=ESTADO_CHOICES, default='configuracion')
    
    # Control de tiempo
    tiempo_limite = models.PositiveIntegerField(default=30, help_text='Tiempo límite en minutos')
    fecha_inicio = models.DateTimeField(null=True, blank=True, help_text='Fecha y hora de inicio programada')
    fecha_inicio_real = models.DateTimeField(null=True, blank=True, help_text='Fecha y hora real de inicio')
    fecha_finalizacion = models.DateTimeField(null=True, blank=True, help_text='Fecha y hora de finalización')
    
    # Configuración avanzada
    mostrar_resultados_inmediatos = models.BooleanField(default=True, help_text='Mostrar resultados al finalizar')
    permitir_reingreso = models.BooleanField(default=False, help_text='Permitir reingresar si se desconecta')
    orden_preguntas_aleatorio = models.BooleanField(default=False, help_text='Mostrar preguntas en orden aleatorio')
    
    # Pregunta actual (para controlar el flujo)
    pregunta_actual = models.ForeignKey('PreguntaCuestionario', null=True, blank=True, 
                                       on_delete=models.SET_NULL, related_name='competencias_actuales')
    
    # Metadatos
    fecha_creacion = models.DateTimeField(auto_now_add=True)
    fecha_modificacion = models.DateTimeField(auto_now=True)
    
    # ============== CLASE META CORREGIDA ==============
    class Meta:
        ordering = ['-fecha_creacion']
        indexes = [
            models.Index(fields=['estado']),
            models.Index(fields=['pin_acceso']),
            models.Index(fields=['codigo_acceso']),
            models.Index(fields=['fecha_inicio_real']),
        ]
    
    # ============== MÉTODOS EXISTENTES ==============
    def save(self, *args, **kwargs):
        if not self.codigo_acceso:
            self.codigo_acceso = generar_codigo_aleatorio()
        if not self.pin_acceso:
            self.pin_acceso = generar_pin_competencia()
        super().save(*args, **kwargs)
    
    def __str__(self):
        return f"Competencia - {self.recurso.titulo}"
    
    def get_url_acceso_estudiante(self):
        """Retorna la URL completa para que estudiantes accedan"""
        from django.urls import reverse
        return reverse('acceder_competencia', kwargs={'pin': self.pin_acceso})
    
    # ============== MÉTODOS CORREGIDOS ==============
    def puede_iniciar(self):
        """Verifica si la competencia puede iniciarse - CORREGIDO"""
        from django.contrib.contenttypes.models import ContentType
        competencia_ct = ContentType.objects.get_for_model(Competencia)
        
        tiene_preguntas = PreguntaCuestionario.objects.filter(
            content_type=competencia_ct,
            object_id=self.id
        ).exists()
        
        return (self.estado in ['configuracion', 'esperando'] and 
                tiene_preguntas and
                self.recurso.titulo.strip())
    
    def esta_activa(self):
        """Verifica si la competencia está actualmente activa"""
        return self.estado == 'activa'
    
    def tiempo_restante_segundos(self):
        """Calcula el tiempo restante en segundos"""
        if not self.fecha_inicio_real or self.estado != 'activa':
            return 0
        
        tiempo_transcurrido = (timezone.now() - self.fecha_inicio_real).total_seconds()
        tiempo_limite_segundos = self.tiempo_limite * 60
        return max(0, tiempo_limite_segundos - tiempo_transcurrido)
    
    def calcular_puntaje_total(self):
        """Calcula el puntaje total de todas las preguntas - CORREGIDO"""
        from django.contrib.contenttypes.models import ContentType
        competencia_ct = ContentType.objects.get_for_model(Competencia)
        
        preguntas = PreguntaCuestionario.objects.filter(
            content_type=competencia_ct,
            object_id=self.id
        )
        
        return preguntas.aggregate(
            total=models.Sum('puntaje')
        )['total'] or 0
    
    # ============== MÉTODOS NUEVOS AGREGADOS ==============
    def get_preguntas(self):
        """Obtiene todas las preguntas de esta competencia"""
        from django.contrib.contenttypes.models import ContentType
        competencia_ct = ContentType.objects.get_for_model(Competencia)
        
        return PreguntaCuestionario.objects.filter(
            content_type=competencia_ct,
            object_id=self.id
        ).order_by('orden')

    def count_preguntas(self):
        """Cuenta el número de preguntas"""
        return self.get_preguntas().count()

    def get_participantes_activos(self):
        """Obtiene participantes activos"""
        return self.participaciones.filter(activo=True)

    def get_participantes_finalizados(self):
        """Obtiene participantes que han finalizado"""
        return self.participaciones.filter(finalizo=True)

    def puede_abrir_sala_espera(self):
        """Verifica si se puede abrir la sala de espera"""
        return (self.estado == 'configuracion' and 
                self.count_preguntas() > 0 and
                self.recurso.titulo.strip())

    def puede_finalizar(self):
        """Verifica si se puede finalizar la competencia"""
        return self.estado in ['esperando', 'activa']

    def esta_en_configuracion(self):
        """Verifica si está en configuración"""
        return self.estado == 'configuracion'

    def esta_esperando(self):
        """Verifica si está esperando participantes"""
        return self.estado == 'esperando'

    def esta_finalizada(self):
        """Verifica si ya finalizó"""
        return self.estado == 'finalizada'

    def get_siguiente_pregunta(self):
        """Obtiene la siguiente pregunta en orden"""
        if not self.pregunta_actual:
            return self.get_preguntas().first()
        
        return self.get_preguntas().filter(
            orden__gt=self.pregunta_actual.orden
        ).first()

    def get_pregunta_anterior(self):
        """Obtiene la pregunta anterior en orden"""
        if not self.pregunta_actual:
            return None
        
        return self.get_preguntas().filter(
            orden__lt=self.pregunta_actual.orden
        ).last()

class GrupoCompetencia(models.Model):
    """Modelo para manejar grupos en competencias grupales"""
    competencia = models.ForeignKey('Competencia', on_delete=models.CASCADE, related_name='grupos')
    nombre = models.CharField(max_length=100, help_text='Nombre del grupo')
    codigo_grupo = models.CharField(max_length=8, help_text='Código único del grupo')
    
    # Control de miembros
    max_miembros = models.PositiveIntegerField(default=4)
    es_completo = models.BooleanField(default=False)
    
    # Metadatos
    fecha_creacion = models.DateTimeField(auto_now_add=True)
    creado_automaticamente = models.BooleanField(default=False)
    
    # ============== CLASE META CORREGIDA ==============
    class Meta:
        unique_together = ['competencia', 'codigo_grupo']
        indexes = [
            models.Index(fields=['competencia']),
            models.Index(fields=['creado_automaticamente']),
        ]
    
    # ============== MÉTODOS EXISTENTES ==============
    def __str__(self):
        return f"Grupo {self.nombre} - {self.competencia.recurso.titulo}"
    
    def save(self, *args, **kwargs):
        if not self.codigo_grupo:
            self.codigo_grupo = ''.join(random.choices(string.ascii_uppercase + string.digits, k=6))
        super().save(*args, **kwargs)
    
    def miembros_count(self):
        """Cuenta los miembros actuales del grupo"""
        return self.participaciones.count()
    
    def puede_agregar_miembro(self):
        """Verifica si se puede agregar un miembro más"""
        return self.miembros_count() < self.max_miembros
    
    def esta_completo(self):
        """Verifica si el grupo está completo"""
        return self.miembros_count() >= self.max_miembros
    
    # ============== MÉTODOS NUEVOS AGREGADOS ==============
    def get_participantes(self):
        """Obtiene todos los participantes del grupo"""
        return self.participaciones.all().select_related('estudiante')

    def get_participantes_activos(self):
        """Obtiene participantes activos del grupo"""
        return self.participaciones.filter(activo=True).select_related('estudiante')

    def calcular_puntaje_promedio(self):
        """Calcula el puntaje promedio del grupo"""
        participaciones = self.participaciones.filter(finalizo=True)
        if not participaciones.exists():
            return 0
        
        from django.db.models import Avg
        return participaciones.aggregate(
            promedio=Avg('puntaje_total')
        )['promedio'] or 0

    def calcular_tiempo_promedio(self):
        """Calcula el tiempo promedio del grupo"""
        participaciones = self.participaciones.filter(
            finalizo=True, 
            tiempo_total_segundos__isnull=False
        )
        if not participaciones.exists():
            return 0
        
        from django.db.models import Avg
        return participaciones.aggregate(
            promedio=Avg('tiempo_total_segundos')
        )['promedio'] or 0

    def todos_finalizaron(self):
        """Verifica si todos los miembros del grupo han finalizado"""
        total_miembros = self.participaciones.filter(activo=True).count()
        finalizados = self.participaciones.filter(finalizo=True).count()
        return total_miembros > 0 and total_miembros == finalizados



class ParticipacionCompetencia(models.Model):
    """Modelo para registrar la participación de estudiantes en competencias"""
    competencia = models.ForeignKey('Competencia', on_delete=models.CASCADE, related_name='participaciones')
    estudiante = models.ForeignKey('User', on_delete=models.CASCADE, related_name='participaciones_competencias')
    
    # Para competencias grupales
    grupo = models.ForeignKey('GrupoCompetencia', null=True, blank=True, 
                             on_delete=models.CASCADE, related_name='participaciones')
    
    # Control de sesión
    fecha_ingreso = models.DateTimeField(auto_now_add=True)
    fecha_ultima_actividad = models.DateTimeField(auto_now=True)
    activo = models.BooleanField(default=True)
    
    # Progreso
    pregunta_actual_orden = models.PositiveIntegerField(default=1)
    finalizo = models.BooleanField(default=False)
    fecha_finalizacion = models.DateTimeField(null=True, blank=True)
    
    # Resultados
    puntaje_total = models.DecimalField(max_digits=8, decimal_places=2, default=0.00)
    tiempo_total_segundos = models.PositiveIntegerField(null=True, blank=True)
    posicion = models.PositiveIntegerField(null=True, blank=True, help_text='Posición en el ranking')
    
    # ============== CLASE META CORREGIDA ==============
    class Meta:
        unique_together = ['competencia', 'estudiante']
        ordering = ['-fecha_ingreso']
        indexes = [
            models.Index(fields=['competencia', 'activo']),
            models.Index(fields=['finalizo']),
            models.Index(fields=['puntaje_total']),
        ]
    
    # ============== MÉTODOS EXISTENTES ==============
    def __str__(self):
        grupo_info = f" (Grupo: {self.grupo.nombre})" if self.grupo else ""
        return f"{self.estudiante.username} - {self.competencia.recurso.titulo}{grupo_info}"
    
    def calcular_puntaje_actual(self):
        """Calcula el puntaje actual del participante"""
        puntaje = self.respuestas.aggregate(
            total=models.Sum('puntaje_obtenido')
        )['total'] or 0
        self.puntaje_total = puntaje
        self.save()
        return puntaje
    
    # ============== MÉTODOS NUEVOS AGREGADOS ==============
    def get_respuestas_correctas(self):
        """Obtiene el número de respuestas correctas"""
        return self.respuestas.filter(es_correcta=True).count()

    def get_respuestas_incorrectas(self):
        """Obtiene el número de respuestas incorrectas"""
        return self.respuestas.filter(es_correcta=False).count()

    def get_respuestas_pendientes(self):
        """Obtiene el número de respuestas pendientes de calificar"""
        return self.respuestas.filter(es_correcta__isnull=True).count()

    def calcular_porcentaje_aciertos(self):
        """Calcula el porcentaje de aciertos"""
        total_respuestas = self.respuestas.count()
        if total_respuestas == 0:
            return 0
        
        correctas = self.get_respuestas_correctas()
        return round((correctas / total_respuestas) * 100, 2)

    def get_tiempo_promedio_respuesta(self):
        """Calcula el tiempo promedio por respuesta"""
        respuestas_con_tiempo = self.respuestas.filter(
            tiempo_respuesta_segundos__isnull=False
        )
        if not respuestas_con_tiempo.exists():
            return 0
        
        from django.db.models import Avg
        return respuestas_con_tiempo.aggregate(
            promedio=Avg('tiempo_respuesta_segundos')
        )['promedio'] or 0

    def esta_en_grupo(self):
        """Verifica si está en un grupo"""
        return self.grupo is not None

    def puede_unirse_grupo(self):
        """Verifica si puede unirse a un grupo"""
        return (self.competencia.modalidad == 'grupal' and 
                self.grupo is None and 
                self.competencia.grupos_abiertos and
                self.competencia.estado in ['configuracion', 'esperando'])


class RespuestaCompetencia(models.Model):
    """Modelo para las respuestas en competencias"""
    participacion = models.ForeignKey('ParticipacionCompetencia', on_delete=models.CASCADE, related_name='respuestas')
    pregunta = models.ForeignKey('PreguntaCuestionario', on_delete=models.CASCADE)
    
    # Tipos de respuesta (similar a RespuestaEstudiante pero simplificado para competencias)
    opcion_seleccionada = models.ForeignKey('OpcionPregunta', on_delete=models.SET_NULL, null=True, blank=True)
    opciones_multiples = models.TextField(null=True, blank=True)
    respuesta_texto = models.TextField(blank=True, null=True)
    respuestas_completar = models.TextField(null=True, blank=True)
    conexiones_realizadas = models.TextField(null=True, blank=True)
    respuesta_smiles = models.TextField(blank=True, null=True)
    
    # Evaluación
    es_correcta = models.BooleanField(null=True, blank=True)
    puntaje_obtenido = models.DecimalField(max_digits=5, decimal_places=2, default=0.00)
    
    # Metadatos de tiempo
    fecha_respuesta = models.DateTimeField(auto_now_add=True)
    tiempo_respuesta_segundos = models.PositiveIntegerField(null=True, blank=True)
    
    class Meta:
        unique_together = ['participacion', 'pregunta']
    
    def __str__(self):
        return f"Respuesta de {self.participacion.estudiante.username} - Pregunta {self.pregunta.id}"
    
    def calcular_puntaje_automatico(self):
        """Calcula el puntaje automáticamente según el tipo de pregunta"""
        tipo = self.pregunta.tipo.nombre
        puntaje_pregunta = self.pregunta.puntaje
        
        if tipo == 'opcion_unica' or tipo == 'falso_verdadero':
            if self.opcion_seleccionada and self.opcion_seleccionada.es_correcta:
                self.puntaje_obtenido = puntaje_pregunta
                self.es_correcta = True
            else:
                self.puntaje_obtenido = 0
                self.es_correcta = False
        
        elif tipo == 'opcion_multiple':
            if self.opciones_multiples:
                ids_seleccionadas = set(self.opciones_multiples.split(','))
                opciones_correctas = set(str(op.id) for op in self.pregunta.opciones.filter(es_correcta=True))
                
                if ids_seleccionadas == opciones_correctas:
                    self.puntaje_obtenido = puntaje_pregunta
                    self.es_correcta = True
                else:
                    self.puntaje_obtenido = 0
                    self.es_correcta = False
        
        # Más tipos se pueden agregar aquí...
        
        self.save()
        return self.puntaje_obtenido

class TipoPregunta(models.Model):
    nombre = models.CharField(max_length=50, unique=True)  # Ej: "seleccion_multiple", "falso_verdadero", "respuesta_abierta", "armar_molecula"
    descripcion = models.TextField(blank=True)

    def __str__(self):
        return self.nombre


class PreguntaCuestionario(models.Model):
    # Relación genérica para que funcione con cuestionarios Y competencias
    content_type = models.ForeignKey(ContentType, on_delete=models.CASCADE)
    object_id = models.PositiveIntegerField()
    content_object = GenericForeignKey('content_type', 'object_id')
    enunciado = models.TextField()
    tipo = models.ForeignKey('TipoPregunta', on_delete=models.PROTECT)
    orden = models.PositiveIntegerField(default=1)
    puntaje = models.PositiveIntegerField(default=1)
    imagen = models.ImageField(upload_to='preguntas/', null=True, blank=True)    
    
    # NUEVOS CAMPOS PARA RECURSOS MULTIMEDIA:
    archivo_multimedia = models.FileField(
        upload_to='preguntas/multimedia/', 
        null=True, 
        blank=True,
        help_text='Archivo multimedia de apoyo (imagen, video, documento)'
    )
    url_youtube = models.CharField(
        max_length=11, 
        null=True, 
        blank=True,
        help_text='ID del video de YouTube (solo el ID, no la URL completa)'
    )
    calificacion_manual = models.BooleanField(
        default=False,
        help_text='Requiere calificación manual'
    )
    fecha_creacion = models.DateTimeField(auto_now_add=True)
    fecha_modificacion = models.DateTimeField(auto_now=True)
    
    # Campos para respuesta abierta
    respuesta_modelo = models.TextField(null=True, blank=True, help_text='Respuesta modelo para guiar evaluación')
    criterios_evaluacion = models.TextField(null=True, blank=True, help_text='Criterios de evaluación')
    longitud_minima = models.PositiveIntegerField(null=True, blank=True, default=50, help_text='Longitud mínima de respuesta')
    longitud_maxima = models.PositiveIntegerField(null=True, blank=True, default=1000, help_text='Longitud máxima de respuesta')
    
    # Campos para completar texto
    texto_completar = models.TextField(null=True, blank=True, help_text='Texto con espacios en blanco marcados con _____')
    respuestas_completar = models.TextField(null=True, blank=True, help_text='Respuestas correctas separadas por comas')
    sensible_mayusculas = models.BooleanField(default=False, help_text='Sensible a mayúsculas y minúsculas')
    ignorar_espacios = models.BooleanField(default=True, help_text='Ignorar espacios adicionales')
    permitir_alternativas = models.BooleanField(default=False, help_text='Permitir respuestas alternativas')
    respuestas_alternativas = models.TextField(null=True, blank=True, help_text='Respuestas alternativas por línea')
    
    # Campos para unir líneas
    columna_izquierda = models.TextField(null=True, blank=True, help_text='Elementos de la columna izquierda, uno por línea')
    columna_derecha = models.TextField(null=True, blank=True, help_text='Elementos de la columna derecha, uno por línea')
    conexiones_correctas = models.TextField(null=True, blank=True, help_text='Conexiones correctas formato: 1-A, 2-B')
    mezclar_opciones = models.BooleanField(default=True, help_text='Mezclar orden de las opciones')
    permitir_conexiones_multiples = models.BooleanField(default=False, help_text='Permitir múltiples conexiones')
    
    # Campos para simuladores
    smiles_objetivo = models.TextField(null=True, blank=True, help_text='SMILES de la molécula objetivo')
    descripcion_molecula = models.TextField(null=True, blank=True, help_text='Descripción de la molécula')
    tolerancia_similitud = models.PositiveIntegerField(null=True, blank=True, default=95, help_text='Porcentaje de similitud requerido')
    permitir_isomeros = models.BooleanField(default=False, help_text='Permitir isómeros estructurales')
    def __str__(self):
        return f"{self.orden}. {self.enunciado[:50]}"

    def get_url_youtube_embed(self):
        """Retorna la URL completa para embed de YouTube"""
        if self.url_youtube:
            return f"https://www.youtube.com/embed/{self.url_youtube}"
        return None
    
    @property
    def cuestionario(self):
        """Retorna el cuestionario si está relacionado con uno"""
        if isinstance(self.content_object, Cuestionario):
            return self.content_object
        return None
    
    @property
    def competencia(self):
        """Retorna la competencia si está relacionada con una"""
        if isinstance(self.content_object, Competencia):
            return self.content_object
        return None
    

class OpcionPregunta(models.Model):
    pregunta = models.ForeignKey('PreguntaCuestionario', on_delete=models.CASCADE, related_name='opciones')
    texto = models.CharField(max_length=255)
    es_correcta = models.BooleanField(default=False)

    def __str__(self):
        return self.texto

class MoleculaEsperadaPregunta(models.Model):
    pregunta = models.OneToOneField('PreguntaCuestionario', on_delete=models.CASCADE, related_name='molecula_objetivo')
    smiles_correcto = models.TextField()
    descripcion_esperada = models.TextField(blank=True)

    def __str__(self):
        return f"SMILES esperado para: {self.pregunta.enunciado[:40]}"


class RespuestaCuestionario(models.Model):
    estudiante = models.ForeignKey('User', on_delete=models.CASCADE, related_name='respuestas_cuestionarios')
    tiempo_respuesta = models.IntegerField(null=True, blank=True)  # En segundos
    pregunta = models.ForeignKey('PreguntaCuestionario', on_delete=models.CASCADE, related_name='respuestas')
    opcion_seleccionada = models.ForeignKey('OpcionPregunta', on_delete=models.SET_NULL, null=True, blank=True)
    respuesta_abierta = models.TextField(blank=True, null=True)
    smiles_respuesta = models.TextField(blank=True, null=True)  # Para preguntas tipo 'armar_molecula'
    fecha_respuesta = models.DateTimeField(auto_now_add=True)
    puntos_obtenidos = models.PositiveIntegerField(default=0)

class IntentoCuestionario(models.Model):
    estudiante = models.ForeignKey('User', on_delete=models.CASCADE, related_name='intentos_cuestionarios')
    cuestionario = models.ForeignKey('Cuestionario', on_delete=models.CASCADE, related_name='intentos')
    numero_intento = models.PositiveIntegerField(default=1)
    fecha_inicio = models.DateTimeField(auto_now_add=True)
    fecha_finalizacion = models.DateTimeField(null=True, blank=True)
    puntaje_obtenido = models.DecimalField(max_digits=6, decimal_places=2, default=0.00)
    tiempo_empleado = models.PositiveIntegerField(null=True, blank=True, help_text='Tiempo en segundos')
    completado = models.BooleanField(default=False)
    
    class Meta:
        unique_together = ['estudiante', 'cuestionario', 'numero_intento']
        ordering = ['-fecha_inicio']
    
    def __str__(self):
        return f"{self.estudiante.username} - {self.cuestionario.recurso.titulo} (Intento {self.numero_intento})"


class ProgresoUsuario(models.Model):
    usuario = models.ForeignKey('User', on_delete=models.CASCADE)
    recurso = models.ForeignKey('Recurso', on_delete=models.CASCADE)
    visto = models.BooleanField(default=False)
    fecha_visto = models.DateTimeField(auto_now_add=True)

    class Meta:
        unique_together = ('usuario', 'recurso')


class RespuestaEstudiante(models.Model):
    intento = models.ForeignKey('IntentoCuestionario', on_delete=models.CASCADE, related_name='respuestas')
    pregunta = models.ForeignKey('PreguntaCuestionario', on_delete=models.CASCADE)
    
    # === RESPUESTAS PARA PREGUNTAS DE OPCIÓN ===
    opcion_seleccionada = models.ForeignKey('OpcionPregunta', on_delete=models.SET_NULL, null=True, blank=True)
    opciones_multiples = models.TextField(null=True, blank=True, help_text='IDs de opciones separadas por comas para opción múltiple')
    
    # === RESPUESTAS PARA PREGUNTAS ABIERTAS ===
    respuesta_texto = models.TextField(blank=True, null=True, help_text='Respuesta de texto libre')
    
    # === RESPUESTAS PARA COMPLETAR TEXTO ===
    respuestas_completar = models.TextField(null=True, blank=True, help_text='Respuestas para completar separadas por |')
    
    # === RESPUESTAS PARA UNIR LÍNEAS ===
    conexiones_realizadas = models.TextField(null=True, blank=True, help_text='Conexiones realizadas formato JSON: {"1":"A", "2":"B"}')
    
    # === RESPUESTAS PARA SIMULADORES ===
    respuesta_smiles = models.TextField(blank=True, null=True, help_text='SMILES para simuladores')
    archivo_molecula = models.FileField(upload_to='respuestas/moleculas/', null=True, blank=True, help_text='Archivo de molécula (.mol, .sdf)')
    descripcion_molecula = models.TextField(null=True, blank=True, help_text='Descripción de la molécula creada')
    
    # === CAMPOS DE EVALUACIÓN ===
    es_correcta = models.BooleanField(null=True, blank=True)
    puntaje_obtenido = models.DecimalField(max_digits=5, decimal_places=2, default=0.00)
    retroalimentacion = models.TextField(null=True, blank=True, help_text='Comentarios del profesor para calificación manual')
    
    # === METADATOS ===
    fecha_respuesta = models.DateTimeField(auto_now_add=True)
    tiempo_respuesta = models.PositiveIntegerField(null=True, blank=True, help_text='Tiempo en segundos para responder esta pregunta')
    
    class Meta:
        unique_together = ['intento', 'pregunta']
        ordering = ['pregunta__orden']
    
    def __str__(self):
        return f"Respuesta de {self.intento.estudiante.username} - Pregunta {self.pregunta.id}"
    
    def get_respuesta_display(self):
        """Retorna una representación legible de la respuesta según el tipo"""
        tipo = self.pregunta.tipo.nombre
        
        if tipo == 'opcion_unica':
            return self.opcion_seleccionada.texto if self.opcion_seleccionada else "Sin respuesta"
        
        elif tipo == 'opcion_multiple':
            if self.opciones_multiples:
                ids = self.opciones_multiples.split(',')
                opciones = OpcionPregunta.objects.filter(id__in=ids)
                return ", ".join([op.texto for op in opciones])
            return "Sin respuesta"
        
        elif tipo == 'falso_verdadero':
            return self.opcion_seleccionada.texto if self.opcion_seleccionada else "Sin respuesta"
        
        elif tipo == 'respuesta_abierta':
            return self.respuesta_texto or "Sin respuesta"
        
        elif tipo == 'completar':
            return self.respuestas_completar or "Sin respuesta"
        
        elif tipo == 'unir_lineas':
            if self.conexiones_realizadas:
                try:
                    import json
                    conexiones = json.loads(self.conexiones_realizadas)
                    return f"Conexiones: {conexiones}"
                except:
                    return self.conexiones_realizadas
            return "Sin conexiones"
        
        elif tipo in ['simulador_2d', 'simulador_3d']:
            return self.respuesta_smiles or "Sin molécula"
        
        return "Respuesta no definida"
    
    def calcular_puntaje_automatico(self):
        """Calcula el puntaje automáticamente según el tipo de pregunta"""
        tipo = self.pregunta.tipo.nombre
        puntaje_pregunta = self.pregunta.puntaje
        
        if tipo == 'opcion_unica':
            if self.opcion_seleccionada and self.opcion_seleccionada.es_correcta:
                self.puntaje_obtenido = puntaje_pregunta
                self.es_correcta = True
            else:
                self.puntaje_obtenido = 0
                self.es_correcta = False
        
        elif tipo == 'opcion_multiple':
            if self.opciones_multiples:
                # Lógica para opción múltiple (todas correctas = puntaje completo)
                ids_seleccionadas = set(self.opciones_multiples.split(','))
                opciones_correctas = set(str(op.id) for op in self.pregunta.opciones.filter(es_correcta=True))
                opciones_incorrectas = set(str(op.id) for op in self.pregunta.opciones.filter(es_correcta=False))
                
                # Solo correctas seleccionadas, ninguna incorrecta
                if ids_seleccionadas == opciones_correctas:
                    self.puntaje_obtenido = puntaje_pregunta
                    self.es_correcta = True
                elif ids_seleccionadas.intersection(opciones_incorrectas):
                    # Seleccionó alguna incorrecta
                    self.puntaje_obtenido = 0
                    self.es_correcta = False
                else:
                    # Parcialmente correcta
                    porcentaje = len(ids_seleccionadas.intersection(opciones_correctas)) / len(opciones_correctas)
                    self.puntaje_obtenido = puntaje_pregunta * porcentaje
                    self.es_correcta = porcentaje >= 1.0
        
        elif tipo == 'falso_verdadero':
            if self.opcion_seleccionada and self.opcion_seleccionada.es_correcta:
                self.puntaje_obtenido = puntaje_pregunta
                self.es_correcta = True
            else:
                self.puntaje_obtenido = 0
                self.es_correcta = False
        
        elif tipo == 'completar':
            if self.respuestas_completar and self.pregunta.respuestas_completar:
                # Comparar respuestas (considerando configuraciones)
                respuestas_correctas = self.pregunta.respuestas_completar.split(',')
                respuestas_estudiante = self.respuestas_completar.split('|')
                
                correctas = 0
                for i, resp_correcta in enumerate(respuestas_correctas):
                    if i < len(respuestas_estudiante):
                        resp_est = respuestas_estudiante[i].strip()
                        resp_cor = resp_correcta.strip()
                        
                        if not self.pregunta.sensible_mayusculas:
                            resp_est = resp_est.lower()
                            resp_cor = resp_cor.lower()
                        
                        if self.pregunta.ignorar_espacios:
                            resp_est = resp_est.replace(' ', '')
                            resp_cor = resp_cor.replace(' ', '')
                        
                        if resp_est == resp_cor:
                            correctas += 1
                
                porcentaje = correctas / len(respuestas_correctas) if respuestas_correctas else 0
                self.puntaje_obtenido = puntaje_pregunta * porcentaje
                self.es_correcta = porcentaje >= 1.0
        
        elif tipo == 'unir_lineas':
            if self.conexiones_realizadas and self.pregunta.conexiones_correctas:
                try:
                    import json
                    conexiones_estudiante = json.loads(self.conexiones_realizadas)
                    conexiones_correctas = {}
                    
                    # Parsear conexiones correctas "1-A, 2-B"
                    for conexion in self.pregunta.conexiones_correctas.split(','):
                        if '-' in conexion:
                            izq, der = conexion.strip().split('-')
                            conexiones_correctas[izq.strip()] = der.strip()
                    
                    correctas = 0
                    for izq, der in conexiones_correctas.items():
                        if conexiones_estudiante.get(izq) == der:
                            correctas += 1
                    
                    porcentaje = correctas / len(conexiones_correctas) if conexiones_correctas else 0
                    self.puntaje_obtenido = puntaje_pregunta * porcentaje
                    self.es_correcta = porcentaje >= 1.0
                except:
                    self.puntaje_obtenido = 0
                    self.es_correcta = False
        
        elif tipo in ['simulador_2d', 'simulador_3d']:
            # Requiere validación manual o comparación de SMILES
            if self.respuesta_smiles and self.pregunta.smiles_objetivo:
                # Comparación básica de SMILES (se puede mejorar con RDKit)
                if self.respuesta_smiles.strip() == self.pregunta.smiles_objetivo.strip():
                    self.puntaje_obtenido = puntaje_pregunta
                    self.es_correcta = True
                else:
                    self.puntaje_obtenido = 0
                    self.es_correcta = False
            else:
                # Requiere calificación manual
                self.es_correcta = None
                self.puntaje_obtenido = 0
        
        elif tipo == 'respuesta_abierta':
            # Siempre requiere calificación manual
            self.es_correcta = None
            self.puntaje_obtenido = 0
        
        self.save()
        return self.puntaje_obtenido

#EVENTOS COMPETENCIAS

class EventoCompetencia(models.Model):
    """Modelo para registrar eventos importantes de la competencia"""
    TIPO_EVENTO_CHOICES = [
        ('creada', 'Competencia Creada'),
        ('configurada', 'Configuración Actualizada'),
        ('sala_abierta', 'Sala de Espera Abierta'),
        ('iniciada', 'Competencia Iniciada'),
        ('pregunta_cambiada', 'Pregunta Cambiada'),
        ('pausada', 'Competencia Pausada'),
        ('reanudada', 'Competencia Reanudada'),
        ('finalizada', 'Competencia Finalizada'),
        ('cancelada', 'Competencia Cancelada'),
        ('estudiante_unido', 'Estudiante se Unió'),
        ('estudiante_salio', 'Estudiante Salió'),
        ('grupo_creado', 'Grupo Creado'),
        ('grupo_eliminado', 'Grupo Eliminado'),
    ]
    
    competencia = models.ForeignKey('Competencia', on_delete=models.CASCADE, related_name='eventos')
    tipo = models.CharField(max_length=20, choices=TIPO_EVENTO_CHOICES)
    descripcion = models.TextField(blank=True)
    usuario = models.ForeignKey('User', on_delete=models.SET_NULL, null=True, blank=True)
    datos_adicionales = models.JSONField(null=True, blank=True)  # Para almacenar datos extra
    fecha_evento = models.DateTimeField(auto_now_add=True)
    
    class Meta:
        ordering = ['-fecha_evento']
    
    def __str__(self):
        return f"{self.competencia.recurso.titulo} - {self.get_tipo_display()} - {self.fecha_evento}"



