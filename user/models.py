from django.db import models
from django.contrib.auth.models import AbstractUser


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

class Curso(models.Model):
    titulo = models.CharField(max_length=200)
    descripcion = models.TextField()
    imagen_portada = models.ImageField(upload_to='cursos/', blank=True, null=True)
    profesor = models.ForeignKey('User', on_delete=models.CASCADE, related_name='cursos_creados')
    categoria = models.ForeignKey('CategoriaCurso', on_delete=models.SET_NULL, null=True, blank=True)
    fecha_creacion = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return self.titulo

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


class Competencia(models.Model):
    recurso = models.OneToOneField('Recurso', on_delete=models.CASCADE, related_name='competencia')
    codigo_acceso = models.CharField(max_length=10, unique=True)
    instrucciones = models.TextField(blank=True)
    fecha_creacion = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return f"Competencia - {self.recurso.titulo}"

class TipoPregunta(models.Model):
    nombre = models.CharField(max_length=50, unique=True)  # Ej: "seleccion_multiple", "falso_verdadero", "respuesta_abierta", "armar_molecula"
    descripcion = models.TextField(blank=True)

    def __str__(self):
        return self.nombre


class PreguntaCuestionario(models.Model):
    cuestionario = models.ForeignKey('Cuestionario', on_delete=models.CASCADE, related_name='preguntas')
    enunciado = models.TextField()
    tipo = models.ForeignKey('TipoPregunta', on_delete=models.PROTECT)
    orden = models.PositiveIntegerField(default=1)
    puntaje = models.PositiveIntegerField(default=1)
    imagen = models.ImageField(upload_to='preguntas/', null=True, blank=True)

    def __str__(self):
        return f"{self.orden}. {self.enunciado[:50]}"

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
    pregunta = models.ForeignKey('PreguntaCuestionario', on_delete=models.CASCADE, related_name='respuestas')
    opcion_seleccionada = models.ForeignKey('OpcionPregunta', on_delete=models.SET_NULL, null=True, blank=True)
    respuesta_abierta = models.TextField(blank=True, null=True)
    smiles_respuesta = models.TextField(blank=True, null=True)  # Para preguntas tipo 'armar_molecula'
    fecha_respuesta = models.DateTimeField(auto_now_add=True)
    puntos_obtenidos = models.PositiveIntegerField(default=0)

    def __str__(self):
        return f"Resp. de {self.estudiante.username} a {self.pregunta.id}"

