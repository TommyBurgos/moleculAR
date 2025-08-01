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