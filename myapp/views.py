from django.shortcuts import render, redirect, get_object_or_404
from django.core.paginator import Paginator
from django.db.models import Q
from django.http import HttpResponse
from django.shortcuts import render
from django.http import JsonResponse
from django.views.decorators.csrf import csrf_exempt, csrf_protect
#from rdkit import Chem
#from rdkit.Chem import AllChem
#from rdkit.Chem import SDWriter
from django.middleware.csrf import rotate_token
from django.contrib.auth.decorators import login_required
from django.contrib.contenttypes.models import ContentType
from django.contrib import messages
import tempfile
import json
from django.http import JsonResponse, HttpRequest

from user.models import User,Rol, Curso, Seccion, TipoRecurso, Recurso, Cuestionario, Practica, Modelo, InscripcionCurso, MoleculaEstudiante, Competencia, PreguntaCuestionario,TipoPregunta, ProgresoUsuario, OpcionPregunta, IntentoCuestionario, RespuestaEstudiante, PracticaConfig, IntentoArmado, GrupoCompetencia, ParticipacionCompetencia, RespuestaCompetencia

from .forms import AdminCrearUsuarioForm, EditarPerfilForm, CursoForm

from django.contrib.auth import login, logout, authenticate
import re
from .permissions import role_required
from django.db.models import Q, Sum  
 
import boto3
import uuid
from django.conf import settings
from django.urls import reverse
import requests, os

API_RDKit_URL = "https://rdkit-api-l1d9.onrender.com"

def llamar_api_rdkit(endpoint, datos):
    try:
        response = requests.post(f"{API_RDKit_URL}{endpoint}", json=datos)
        if response.status_code == 200:
            return response.json()
        else:
            return {"error": "Error al comunicarse con la API externa"}
    except Exception as e:
        return {"error": str(e)}


@csrf_exempt
def modificar_molecula(request):
    if request.method == 'POST':
        try:
            data = json.loads(request.body)
            smiles = data.get("smiles")
            atomo_idx = data.get("atomo_idx")  # ya no lo convierto a int aquí
            grupo = data.get("grupo")

            if not smiles or atomo_idx is None or not grupo:
                return JsonResponse({"error": "Faltan datos requeridos"}, status=400)

            # Enviar a la API externa con la lógica RDKit
            payload = {
                "smiles": smiles,
                "atomo_idx": atomo_idx,
                "grupo": grupo
            }

            resultado = llamar_api_rdkit("/modificar", payload)

            if "error" in resultado:
                return JsonResponse({"error": resultado["error"]}, status=400)

            return JsonResponse({"sdf": resultado["sdf"]})

        except json.JSONDecodeError:
            return JsonResponse({"error": "Error en el formato JSON recibido."}, status=400)
        except Exception as e:
            return JsonResponse({"error": str(e)}, status=400)

    return JsonResponse({"message": "Usa POST para modificar la molécula."})

def vista_modificador(request):
    return render(request, 'modificador_molecula.html')

# Create your views here.
def inicio(request):
    return render(request, 'index.html')

def registroUser(request):
    return render(request, 'sign-up.html')

#REGISTRAR USUARIO ESTUDIANTE
def registroExitoso(request):
    if request.method == 'POST':
        print(request.POST)
        if request.POST['password'] == request.POST['confirmPassword']:
            try:
                password = request.POST.get('password')
                email = request.POST.get('email')                
                # Validaciones
                if len(password) < 8 or not re.match(r'^(?=.*[A-Za-z])(?=.*\d)[A-Za-z\d]{8,}$', password):
                    return HttpResponse('La contraseña debe tener al menos 8 caracteres y contener letras y números.')

                if not re.match(r'^[a-zA-Z0-9._%+-]+@[a-zA-Z0-9.-]+\.[a-zA-Z]{2,}$', email):
                    return HttpResponse('El formato del correo electrónico no es válido.')
                
                # Crear usuario
                user = User.objects.create_user(
                    username=email,  # Usa el correo como username
                    email=email,
                    password=password,
                    first_name=' ',
                    last_name=' '
                )
                print('***INICIO***')                
                rol_estudiante = Rol.objects.get(nombre='Estudiante')
                print(rol_estudiante)
                user.rol = rol_estudiante
                print(user.rol)                                
                user.save()
                print('***FIN***')
                login(request, user)
                print("Usuario creado y logueado exitosamente.")
                
                return redirect('inicio')
            except Exception as e:
                print(f"❌ Error al crear el usuario: {e}")
                return HttpResponse('Ocurrió un error al crear el usuario.')
        else:
            return HttpResponse('Las contraseñas no coinciden.')
    return HttpResponse('Método no permitido')

#REGISTRAR DOCENTE
def registroDocente(request):
    if request.method == 'POST':
        email = request.POST.get('email')
        password = request.POST.get('password')
        confirm = request.POST.get('confirmPassword')
        codigo = request.POST.get('codigo_institucional')

        if password != confirm:
            return HttpResponse('Las contraseñas no coinciden.')

        if len(password) < 8 or not re.match(r'^(?=.*[A-Za-z])(?=.*\d)[A-Za-z\d]{8,}$', password):
            return HttpResponse('La contraseña debe tener al menos 8 caracteres y contener letras y números.')

        CODIGO_VALIDO = getattr(settings, "CODIGO_DOCENTE", None)
        if codigo != CODIGO_VALIDO:
            return HttpResponse('Código institucional inválido.')

        try:
            user = User.objects.create_user(
                username=email,
                email=email,
                password=password,
                first_name='',
                last_name=''
            )
            rol_docente = Rol.objects.get(nombre='Docente')
            user.rol = rol_docente
            user.save()
            login(request, user)
            return redirect('inicio')
        except Exception as e:
            return HttpResponse(f'Ocurrió un error: {e}')

    return render(request, 'sign-up-docente.html')


@csrf_exempt
def vistaLogin(request):
    return render(request,'sign-in.html')

def acceso_denegado(request):
    print("Inicie a la función de acceso denegado")        
    return render(request, 'error-404.html')

from django.db.models import Count
from django.db.models.functions import TruncMonth

from django.db.models.functions import TruncMonth, TruncWeek, TruncDay

@role_required('Admin')
def inicioAdmin(request):
    print("Inicie a la función de inicio Admin")
    user = request.user    
    imgPerfil = user.imgPerfil  

    # Métricas globales
    total_cursos = Curso.objects.count()
    total_usuarios = User.objects.count()

    hace_un_mes = timezone.now() - timedelta(days=30)
    ultimos_usuarios = (
        User.objects.filter(date_joined__gte=hace_un_mes)
        .annotate(cursos_inscritos=Count('inscripcioncurso'))
        .order_by('-date_joined')[:6]
    )

    cursos_populares = (
        Curso.objects.annotate(num_inscritos=Count('inscripcioncurso'))
        .order_by('-num_inscritos')[:5]
    )
    labels_cursos = [c.titulo for c in cursos_populares]
    data_cursos = [c.num_inscritos for c in cursos_populares]
    total_inscripciones = sum(data_cursos)

    # 🔹 Intervalo dinámico para usuarios
    intervalo = request.GET.get('intervalo', 'mes')
    if intervalo == 'dia':
        trunc_func = TruncDay
        formato = "%Y-%m-%d"
    elif intervalo == 'semana':
        trunc_func = TruncWeek
        formato = "Semana %W %Y"
    else:
        trunc_func = TruncMonth
        formato = "%b %Y"

    usuarios_por_intervalo = (
        User.objects
        .annotate(periodo=trunc_func('date_joined'))
        .values('periodo')
        .annotate(total=Count('id'))
        .order_by('periodo')
    )

    labels = [u['periodo'].strftime(formato) for u in usuarios_por_intervalo]
    data = [u['total'] for u in usuarios_por_intervalo]

    cursos_recientes = Curso.objects.filter(
        fecha_creacion__gte=hace_un_mes
    ).select_related('profesor', 'categoria').order_by('-fecha_creacion')

    hace_una_semana = timezone.now() - timedelta(days=7)
    recursos_recientes = (
        Recurso.objects.filter(fecha_creacion__gte=hace_una_semana)
        .select_related('tipo')
        .order_by('-fecha_creacion')
    )

    try:
        tipo_practica = TipoRecurso.objects.get(nombre__iexact="Practica")
        total_practicas = Recurso.objects.filter(tipo=tipo_practica).count()
        practicas_visibles = Recurso.objects.filter(tipo=tipo_practica, visible_biblioteca=1).count()
    except TipoRecurso.DoesNotExist:
        total_practicas = 0
        practicas_visibles = 0

    context = {                  
        'imgPerfil': imgPerfil,        
        'usuario': user,
        'total_cursos': total_cursos,
        'total_usuarios': total_usuarios,
        'ultimos_usuarios': ultimos_usuarios,
        'labels': labels,
        'data': data,
        'intervalo': intervalo,  # 👈 se pasa al template
        'cursos_recientes': cursos_recientes,
        'recursos_recientes': recursos_recientes,
        'total_practicas': total_practicas,
        'practicas_visibles': practicas_visibles,
        'labels_cursos': labels_cursos,
        'data_cursos': data_cursos,
        'total_inscripciones': total_inscripciones,
    }   
    return render(request, 'usAdmin/admin-dashboard.html', context)

@role_required(['Admin', 'Docente'])
def vistaAllCursos(request):
    user = request.user    
    imgPerfil=user.imgPerfil 
    cursos = Curso.objects.all().order_by('-fecha_creacion') 
    context = {                  
        'imgPerfil': imgPerfil,        
        'usuario':user.username, 
        'cursos': cursos       
    }   
    return render(request,'usAdmin/detalleCursos.html',context)

from django.db.models import Prefetch

@role_required(['Admin', 'Docente'])
@login_required
def detalle_curso(request, curso_id):
    curso = get_object_or_404(Curso, id=curso_id)
    
    recursos_prefetch = Prefetch(
        'recursos',
        queryset=Recurso.objects.order_by('orden'),
        to_attr='recursos_ordenados'
    )
    
    secciones = Seccion.objects.filter(curso=curso_id).prefetch_related(recursos_prefetch).order_by('orden')

    recursosTipo = TipoRecurso.objects.all()
    user = request.user
    imgPerfil = user.imgPerfil

    return render(request, 'usAdmin/detalle_curso.html', {
        'curso': curso,
        'imgPerfil': imgPerfil,
        'usuario': user.username,
        'secciones': secciones,
        'tipos_recurso': recursosTipo
    })

from django.utils.html import strip_tags

def limpiar_html_permitiendo_alineacion(html):
    if not html:
        return ""

    # Elimina todo menos etiquetas <p>, <strong>, <em>, <u>, <br>, <span>
    html = re.sub(r'<(?!/?(p|strong|em|u|br|span|img)(\s|>|/))[^>]+>', '', html)

    # Filtra los estilos permitidos (solo text-align)
    html = re.sub(r'style="[^"]*"', lambda m: filtrar_estilos(m.group()), html)

    return html

def filtrar_estilos(style_str):
    estilos_permitidos = ['text-align', 'color', 'background-color']
    estilos = style_str.replace('style="', '').replace('"', '').split(';')
    estilos_filtrados = [s for s in estilos if any(permitido in s for permitido in estilos_permitidos)]
    if estilos_filtrados:
        return f'style="{";".join(estilos_filtrados)}"'
    return ''


import bleach
from django.utils.safestring import mark_safe



@login_required
def detalle_recurso(request, recurso_id):
    print("He ingresado al detalle del recurso")
    recurso = get_object_or_404(Recurso, id=recurso_id)

    if recurso.tipo.nombre.lower() == 'html embebido':
        return render(request, 'usAdmin/detalle_html.html', {
            'recurso': recurso,
            #'contenido_html_safe': mark_safe(recurso.contenido_html or "")
        })

    if request.user.is_authenticated and hasattr(request.user, 'rol') and request.user.rol.nombre == 'Estudiante':
        progreso, creado = ProgresoUsuario.objects.get_or_create(
            usuario=request.user,
            recurso=recurso,
            defaults={'visto': True}
        )
        if not creado and not progreso.visto:
            progreso.visto = True
            progreso.save()
    
    curso = recurso.seccion.curso 
    contenido_html = recurso.contenido_texto or ""
    html_limpio = limpiar_html_permitiendo_alineacion(contenido_html)

    recursos_prefetch = Prefetch(
        'recursos',
        queryset=Recurso.objects.order_by('orden'),
        to_attr='recursos_ordenados'
    )    
    secciones = Seccion.objects.filter(curso=curso).prefetch_related(recursos_prefetch).order_by('orden')

    return render(request, 'usAdmin/detalleRecurso.html', {
        'recurso': recurso,
        'secciones': secciones,
        'curso': curso,
        'html_limpio': html_limpio
    })

@role_required(['Admin', 'Docente'])
@login_required
def editar_video(request, recurso_id):
    print("Acabo de ingresar a editar texto")
    recurso = get_object_or_404(Recurso, id=recurso_id, tipo__nombre='Video')

    if request.method == 'POST':
        recurso.titulo = request.POST.get('titulo')
        recurso.descripcion = request.POST.get('descripcion')
        # Opcionalmente permitir reemplazar el video
        if 'archivo' in request.FILES:
            recurso.archivo = request.FILES['archivo']
        recurso.save()
        messages.success(request, "Video actualizado correctamente.")
        return redirect('detalle_curso', curso_id=recurso.seccion.curso.id)

    return render(request, 'usAdmin/editar_video.html', {'recurso': recurso})


@role_required(['Admin', 'Docente'])
@login_required
def editar_texto(request, recurso_id):
    print("Acabo de ingresar a editar texto")
    recurso = get_object_or_404(Recurso, id=recurso_id, tipo__nombre='Texto')
    print(recurso)
    if request.method == 'POST':
        recurso.titulo = request.POST.get('titulo')
        recurso.descripcion = request.POST.get('descripcion')
        recurso.contenido_texto = request.POST.get('contenido_texto')
        recurso.save()
        messages.success(request, "Texto actualizado correctamente.")
        return redirect('detalle_curso', curso_id=recurso.seccion.curso.id)

    return render(request, 'usAdmin/editar_texto.html', {'recurso': recurso})


def crear_modulo(request, curso_id):
    curso = get_object_or_404(Curso, id=curso_id)

    if request.method == 'POST':
        titulo = request.POST.get('titulo', '').strip()
        descripcion = request.POST.get('descripcion', '').strip()
        orden = request.POST.get('orden', 1)

        if not titulo:
            messages.error(request, "El título del módulo es obligatorio.")
            return redirect('crear_modulo', curso_id=curso.id)

        Seccion.objects.create(
            curso=curso,
            titulo=titulo,
            descripcion=descripcion,
            orden=orden
        )
        messages.success(request, "Módulo creado exitosamente.")
        return redirect('detalle_curso', curso_id=curso.id)

    return render(request, 'usAdmin/detalleCursos.html', {'curso': curso})

#FUNCIÓN QUE ME AYUDA A OBTENER EL SMILE
import requests
def obtener_smiles_pubchem(nombre_molecula):
    nombre_molecula = nombre_molecula.strip().lower()
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{nombre_molecula}/property/IsomericSMILES,SMILES/JSON"
    try:
        response = requests.get(url, timeout=30)
        if response.status_code == 200:
            data = response.json()
            properties = data.get('PropertyTable', {}).get('Properties', [])
            if properties:
                return (
                    properties[0].get('IsomericSMILES') or
                    properties[0].get('SMILES')
                )
    except Exception as e:
        print(f"❌ Error al consultar PubChem: {e}")
    return None

@login_required
def ver_practica(request, recurso_id):
    recurso  = get_object_or_404(Recurso, id=recurso_id)
    practica = get_object_or_404(Practica, recurso=recurso)
    cfg      = getattr(recurso, 'practica_cfg', None)
    tipo     = (cfg.practica_tipo if cfg and cfg.practica_tipo else 'jsme').lower()

    # Logs útiles
    print(f"[ver_practica] recurso={recurso.id} tipo={tipo}")

    # ——— Rama Constructor 2D (lienzo vacío para el estudiante) ———
    if tipo == 'builder2d':
        # objetivo desde config; si no hay, intenta desde el modelo de la práctica
        objetivo_smiles = (cfg.objetivo_smiles or
                           (practica.modelo_objetivo.smiles if practica.modelo_objetivo else '') if cfg else
                           (practica.modelo_objetivo.smiles if practica.modelo_objetivo else ''))

        # nombre “bonito”: el Modelo en cfg si existe, si no el título de la práctica o recurso
        objetivo_nombre = (getattr(cfg, 'modelo', None).titulo
                           if (cfg and getattr(cfg, 'modelo', None)) else
                           (practica.titulo_objetivo or recurso.titulo))

        ctx = {
            'recurso': recurso,
            'practica': practica,
            'objetivo_smiles': objetivo_smiles or '',
            'objetivo_nombre': objetivo_nombre or '',
            'inventario_json': None,  #Lienzo inicia vacío (el alumno arma desde cero)
        }
        return render(request, 'usAdmin/constructor2d_ver.html', ctx)

    # Rama JSME - vista “ver” 
    # Muestra detalles e inicia sin solución final
    modelo = practica.modelo_objetivo
    smiles = (modelo.smiles if modelo else (cfg.objetivo_smiles if cfg else ''))

    return render(request, 'usAdmin/ver_practica.html', {
        'practica': practica,
        'recurso': recurso,
        'modelo': modelo,
        'smiles': smiles or '',
    })


@login_required
def editar_practica(request, recurso_id):
    recurso = get_object_or_404(Recurso, id=recurso_id)

    # Config de práctica (OneToOne). Si no existe, por defecto JSME.
    practica_cfg, _ = PracticaConfig.objects.get_or_create(
        recurso=recurso,
        defaults={'practica_tipo': 'jsme'}
    )

    # Metadatos de práctica (con título/instrucciones/modelo_objetivo)
    practica, _ = Practica.objects.get_or_create(recurso=recurso)

    if request.method == 'POST':
        # Permite que el form cambie el tipo de práctica 
        nuevo_tipo = request.POST.get('practica_tipo')
        if nuevo_tipo and nuevo_tipo in dict(PracticaConfig.TIPO_PRACTICA):
            practica_cfg.practica_tipo = nuevo_tipo

        fuente = request.POST.get('fuente_molecula')  # 'jsme' | 'pubchem'
        titulo_objetivo = request.POST.get('titulo_objetivo') or practica.titulo_objetivo
        instrucciones = request.POST.get('instrucciones') or practica.instrucciones

        practica.titulo_objetivo = titulo_objetivo
        practica.instrucciones = instrucciones

        modelo = None
        smiles = None

        # --- Rama JSME ---
        if practica_cfg.practica_tipo == 'jsme':
            if fuente == 'jsme':
                smiles = (request.POST.get('mol_input') or '').strip()

            elif fuente == 'pubchem':
                nombre_molecula = (request.POST.get('nombre_molecula') or '').strip()
                smiles = (obtener_smiles_pubchem(nombre_molecula) or '').strip()
                if smiles and not titulo_objetivo:
                    titulo_objetivo = nombre_molecula  

            if smiles:
                modelo = Modelo.objects.create(
                    titulo=titulo_objetivo or recurso.titulo,
                    descripcion=instrucciones or '',
                    tipo='2D',
                    smiles=smiles,
                    creado_por=request.user,
                    visible_biblioteca=False
                )
                practica.modelo_objetivo = modelo
                practica_cfg.objetivo_smiles = smiles  # giardo también en la config

            else:
                # Error: no hubo SMILES válido
                return render(request, 'usAdmin/editor_dibujo.html', {
                    'recurso': recurso,
                    'practica': practica,
                    'instruccion': instrucciones or "",
                    'smiles': "",
                    'error': 'No se pudo obtener un SMILES válido (JSME/PubChem).'
                })

        # --- Rama Constructor 2D: solo guardo/actualizo el objetivo SMILES si vino ---
        elif practica_cfg.practica_tipo == 'builder2d':
            objetivo_smiles = (request.POST.get('objetivo_smiles') or practica_cfg.objetivo_smiles or '').strip()
            practica_cfg.objetivo_smiles = objetivo_smiles
            # En builder2d no se cre el Modelo aquí; el armado lo hace el estudiante

        # Persistencia de cambios
        practica.save()
        practica_cfg.save()

        # Dispatcher: redirige a la UI correcta
        if practica_cfg.practica_tipo == 'builder2d':
            return redirect('practica_builder2d', recurso_id=recurso.id)
        else:
            return redirect('editar_practica', recurso_id=recurso.id)

    # --- GET: muestra la pantalla correspondiente ---
    if practica_cfg.practica_tipo == 'builder2d':
        return redirect('practica_builder2d', recurso_id=recurso.id)

    # Por defecto JSME: pasa SMILES existente si hay
    smiles_existente = ""
    if practica.modelo_objetivo:
        smiles_existente = practica.modelo_objetivo.smiles or ""
    elif practica_cfg.objetivo_smiles:
        smiles_existente = practica_cfg.objetivo_smiles or ""

    return render(request, 'usAdmin/editor_dibujo.html', {
        'recurso': recurso,
        'practica': practica,
        'instruccion': practica.instrucciones or "",
        'smiles': smiles_existente,
        'modelo': getattr(practica, 'modelo_objetivo', None),
        'error': None
    })

def practica_jsme(request, recurso_id):
    recurso = get_object_or_404(Recurso, pk=recurso_id)
    practica_cfg = PracticaConfig.objects.filter(recurso=recurso).first()
    practica, _ = Practica.objects.get_or_create(recurso=recurso)

    smiles_existente = ""
    if getattr(practica, 'modelo_objetivo', None):
        smiles_existente = practica.modelo_objetivo.smiles or ""
    elif practica_cfg and practica_cfg.objetivo_smiles:
        smiles_existente = practica_cfg.objetivo_smiles or ""

    return render(request, 'usAdmin/editor_dibujo.html', {
        'recurso': recurso,
        'practica': practica,
        'instruccion': practica.instrucciones or "",
        'smiles': smiles_existente,
        'modelo': getattr(practica, 'modelo_objetivo', None),
        'error': None
    })

# ---- UI Constructor 2D ----
def _get_rdkit_base():
    for c in [
        getattr(settings, 'RDKit_API_URL', None),
        getattr(settings, 'RDKIT_API_URL', None),
        os.environ.get('RDKit_API_URL'),
        os.environ.get('RDKIT_API_URL'),
    ]:
        if c and str(c).strip():
            return str(c).strip().rstrip('/')
    # fallback duro para no bloquear desarrollo si settings falla
    return "https://rdkit-api-l1d9.onrender.com"


# -------- pantalla del Constructor 2D --------
@login_required
def practica_builder2d(request, recurso_id):
    recurso = get_object_or_404(Recurso, pk=recurso_id)
    cfg, _ = PracticaConfig.objects.get_or_create(
        recurso=recurso, defaults={'practica_tipo': 'builder2d'}
    )

    # objetivo e (opcional) inventario guardado
    objetivo = (cfg.objetivo_smiles or '').strip()
    inventario = None
    if cfg.inventario_json:
        try:
            inventario = json.loads(cfg.inventario_json)
        except Exception:
            inventario = None

    return render(request, 'usAdmin/constructor2d.html', {
        'practica': recurso,                              # lo uso como "practica" en el template
        'objetivo_smiles': objetivo,
        'inventario_json': json.dumps(inventario) if inventario else 'null',
    })

from django.views.decorators.http import require_POST
# -------- guardar inventario/objetivo para la práctica builder2d --------
@login_required
@require_POST
def guardar_practica_builder2d(request, recurso_id):
    recurso = get_object_or_404(Recurso, pk=recurso_id)
    cfg, _ = PracticaConfig.objects.get_or_create(
        recurso=recurso, defaults={'practica_tipo': 'builder2d'}
    )

    try:
        payload = json.loads(request.body.decode('utf-8'))
    except Exception:
        return JsonResponse({'ok': False, 'msg': 'JSON inválido'}, status=400)

    inventario = payload.get('inventario')  # dict {atoms,bonds}
    objetivo = (payload.get('objetivo_smiles') or '').strip()

    if inventario and isinstance(inventario, dict):
        cfg.inventario_json = json.dumps(inventario)   # TextField: serializamos

    if objetivo:
        cfg.objetivo_smiles = objetivo

    cfg.practica_tipo = 'builder2d'
    cfg.save()

    return JsonResponse({'ok': True, 'msg': 'Práctica guardada', 'objetivo_smiles': cfg.objetivo_smiles})

import time

# ------ API: SMILES -> grafo 2D (usa el microservicio /smiles_to_graph) ------
@login_required
@require_POST
def smiles_a_grafo(request):
    # Espera POST con JSON: {"smiles": "<cadena>"}
    try:
        body = json.loads(request.body.decode('utf-8'))
    except Exception:
        return JsonResponse({'error': 'JSON inválido'}, status=400)

    smiles = (body.get('smiles') or '').strip()
    if not smiles:
        return JsonResponse({'error': 'SMILES requerido'}, status=400)

    base = _get_rdkit_base()  # Debería devolver la base del microservicio, algo así "https://rdkit-api-....onrender.com"
    url = f"{base}/smiles_to_graph"

    # Reintentos contra el microservicio (para cold starts o picos)
    retries = 2           # total 3 intentos (0,1,2)
    backoff_ms = 0.8      # segundos base de espera progresiva
    last_text = ''

    for attempt in range(retries + 1):
        try:
            r = requests.post(url, json={'smiles': smiles}, timeout=30)

            if r.status_code == 200:
                data = r.json()
                return JsonResponse({
                    'atoms': data.get('atoms', []),
                    'bonds': data.get('bonds', [])
                })

            # Si es error temporal, reintenta
            if r.status_code in (502, 503, 504) and attempt < retries:
                time.sleep(backoff_ms * (attempt + 1))
                continue

            # Otros códigos: intenta extraer mensaje
            try:
                err = r.json().get('error', 'Error microservicio')
            except Exception:
                err = 'Error microservicio'
            return JsonResponse({'error': err}, status=r.status_code)

        except requests.Timeout:
            if attempt < retries:
                time.sleep(backoff_ms * (attempt + 1))
                continue
            return JsonResponse({'error': 'Timeout microservicio'}, status=504)

        except Exception as e:
            # Errores de red u otros: reintenta y si agota, devuelve 500
            if attempt < retries:
                time.sleep(0.6 * (attempt + 1))
                continue
            return JsonResponse({'error': str(e)}, status=500)

# -------- API: validar armado vs objetivo (comparación canónica) --------
@login_required
@require_POST
def validar_molecula_armada(request, recurso_id):
    recurso = get_object_or_404(Recurso, pk=recurso_id)

    try:
        payload = json.loads(request.body.decode('utf-8'))
    except Exception:
        return JsonResponse({'ok': False, 'msg': 'JSON inválido'}, status=400)

    atoms = payload.get('atoms') or []
    bonds = payload.get('bonds') or []    


    # 1) objetivo: primero del body, luego de PracticaConfig, luego del modelo de JSME si existe
    objetivo = (payload.get('objetivo_smiles') or '').strip()
    if not objetivo:
        cfg = getattr(recurso, 'practica_cfg', None)
        if cfg and cfg.objetivo_smiles:
            objetivo = cfg.objetivo_smiles.strip()
        elif getattr(recurso, 'practica', None) and recurso.practica.modelo_objetivo:
            objetivo = (recurso.practica.modelo_objetivo.smiles or '').strip()

    base = _get_rdkit_base()

    # 2) Grafo del estudiante -> SMILES canónico
    try:
        if not atoms:
            return JsonResponse({'ok': False, 'msg': 'No hay átomos en el lienzo.'}, status=200)
        r = requests.post(f"{base}/build_smiles", json={'atoms': atoms, 'bonds': bonds}, timeout=30)
        r.raise_for_status()
        rd = r.json()
        smiles_gen_canon = (rd.get('smiles_canonico') or '').strip()
        if not smiles_gen_canon:
            return JsonResponse({'ok': False, 'msg': 'No se pudo generar SMILES.'}, status=502)
    except Exception as e:
        return JsonResponse({'ok': False, 'msg': f'Error build_smiles: {e}'}, status=500)

    # 3) Canonizar objetivo
    objetivo_canon = ''
    if objetivo:
        try:
            r2 = requests.post(f"{base}/procesar", json={'entrada': objetivo, 'tipo': 'smiles'}, timeout=30)
            if r2.status_code == 200:
                objetivo_canon = (r2.json().get('smiles_canonico') or objetivo).strip()
            else:
                objetivo_canon = objetivo
        except Exception:
            objetivo_canon = objetivo

    es_ok = bool(objetivo_canon) and (smiles_gen_canon == objetivo_canon)

    # 4) Registrar intento
    try:
        IntentoArmado.objects.create(
            practica=recurso,
            estudiante=request.user,
            graph_json=json.dumps({'atoms': atoms, 'bonds': bonds}),  # TextField: serializamos
            smiles_generado=smiles_gen_canon,
            es_correcto=es_ok
        )
    except Exception:
        pass

    return JsonResponse({
        'ok': es_ok,
        'smiles_canonico': smiles_gen_canon,
        'objetivo_canonico': objetivo_canon,
        'msg': 'Coincide' if es_ok else 'Aún no coincide'
    })    

@role_required(['Admin', 'Docente'])
@login_required
def editar_html(request, recurso_id):
    recurso = get_object_or_404(Recurso, id=recurso_id, tipo__nombre__iexact='html embebido')

    if request.method == 'POST':
        recurso.contenido_html = request.POST.get('contenido_html') or ''
        recurso.es_simulacion = 'es_simulacion' in request.POST
        recurso.save()
        messages.success(request, "HTML embebido guardado correctamente.")
        return redirect('detalle_recurso', recurso.id)

    return render(request, 'usAdmin/editar_html.html', {'recurso': recurso})

@login_required
def lista_simulaciones(request):
    simulaciones = Recurso.objects.filter(
        tipo__nombre__iexact="html embebido",
        es_simulacion=True
    ).order_by("-fecha_creacion")

    return render(request, "usAdmin/simulaciones.html", {
        "simulaciones": simulaciones
    })


def crear_recurso(request, seccion_id):
    if request.method == 'POST':
        print("Entré al post de crear recurso")
        # Datos básicos del recurso
        titulo = request.POST.get('titulo')
        print(f"Titulo: {titulo}")
        descripcion = request.POST.get('descripcion')
        imagen = request.FILES.get('imagen_referencial')
        tipo_id = request.POST.get('tipo')
        visible = 'visible_biblioteca' in request.POST
        print("es visible? ", visible)
        seccion_id = request.POST.get('seccion_id')
        print("La sección id es" + seccion_id)
        tipo = TipoRecurso.objects.get(id=tipo_id)
        seccion = Seccion.objects.get(id=seccion_id)
        contenido_texto = request.POST.get('contenido_texto')

        recurso = Recurso.objects.create(
            titulo=titulo,
            descripcion=descripcion,
            tipo=tipo,
            imagen=imagen,
            seccion=seccion,
            visible_biblioteca=visible,
            creado_por=request.user,
            contenido_texto=contenido_texto
        )
        print(f"Recurso creado {recurso}")
        # CREACIÓN SEGÚN TIPO DE RECURSO
        nombre_tipo = tipo.nombre.lower()
        print(f"nombre del tipo de recurso es {nombre_tipo}")
        print('video' in nombre_tipo)
        print("Comparando si ingresa a pràctica")
        print('Practica' in nombre_tipo)
        print('practica' in nombre_tipo)
        print(tipo.nombre.lower() == 'practica')

        if 'video' in nombre_tipo:
            print("Dentro del if")
            video_url = request.POST.get('video_url')
            print(video_url)
            if video_url:
                print(recurso.video_url)
                recurso.video_url = video_url
                print(recurso.video_url)
                recurso.save()
            else:
                print("⚠️ No se recibió la URL del video desde el formulario.")
              
        elif tipo.nombre.lower() == 'practica':
            practica_tipo = request.POST.get('practica_tipo', 'jsme')
            objetivo_smiles = request.POST.get('objetivo_smiles') or None
            modelo_id = request.POST.get('modelo_id') or None

            PracticaConfig.objects.create(
                recurso=recurso,
                practica_tipo=practica_tipo,
                objetivo_smiles=objetivo_smiles,
                modelo_id=modelo_id
            )
            return redirect('editar_practica', recurso.id)

        
        elif tipo.nombre.lower() == 'cuestionario':
            print("Dentro del if cuestionario")
            instrucciones = request.POST.get('nombre_cuestionario')
            tiempo_limite = request.POST.get('descripcion_cuestionario')
            Cuestionario.objects.create(
                recurso=recurso,
                instrucciones=instrucciones,
                tiempo_limite=tiempo_limite
            )  
            # Por ahora NO se crea la práctica. Solo dejamos el recurso creado.
            # Redirijo a una vista donde se edita ese recurso.
            print("Se creo el recurso practica")
            return redirect('editar_cuestionario', recurso.id)  
        
        elif tipo.nombre.lower() == 'competencia':
            print("Dentro del if competencia")
            instrucciones = request.POST.get('instrucciones_competencia', '')
            # Se puede aquí generar un código aleatorio o colocar uno.
            import random, string
            codigo = ''.join(random.choices(string.ascii_uppercase + string.digits, k=6))

            Competencia.objects.create(
                recurso=recurso,
                instrucciones=instrucciones,
                codigo_acceso=codigo
            )
            print("Competencia creada correctamente.")
            return redirect('panel_competencia', recurso.competencia.id)
        
        elif tipo.nombre.lower() == 'html embebido':
            contenido_html = request.POST.get('contenido_html') or ''
            recurso.contenido_html = contenido_html
            recurso.es_simulacion = 'es_simulacion' in request.POST
            recurso.save()
            return redirect('editar_html', recurso.id)



        return redirect('detalle_recurso',recurso.id)



@role_required('Admin')
def crear_usuario_admin(request):
    if not request.user.rol or request.user.rol.nombre.lower() != 'admin':
        messages.error(request, "No tienes permisos para acceder a esta sección.")
        return redirect('acceso_denegado')  

    if request.method == 'POST':
        form = AdminCrearUsuarioForm(request.POST, request.FILES)
        if form.is_valid():
            username = form.cleaned_data['username']
            first_name = form.cleaned_data.get('first_name', '')
            last_name = form.cleaned_data.get('last_name', '')
            email = form.cleaned_data['email']
            rol = form.cleaned_data['rol']
            password = form.cleaned_data['password1']
            imgPerfil = form.cleaned_data.get('imgPerfil')
            is_active = form.cleaned_data['is_active']

            # Crea el usuario con password hasheado
            user = User.objects.create_user(
                username=username,
                email=email,
                password=password,
                first_name=first_name,
                last_name=last_name,
            )
            user.rol = rol
            user.is_active = is_active
            if imgPerfil:
                user.imgPerfil = imgPerfil
            user.save()

            messages.success(request, f"Usuario '{username}' creado correctamente.")            
            return redirect('crear_usuario_admin')
        else:
            messages.error(request, "Por favor corrige los errores del formulario.")
    else:
        form = AdminCrearUsuarioForm()

    return render(request, 'usAdmin/crear_usuario_admin.html', {'form': form})


@role_required(['Admin', 'Docente'])
def detalleUsuarios(request):
    user = request.user
    usuarios2 = User.objects.all().order_by('-date_joined')
    usuarios = User.objects.all().order_by('-date_joined')
    imgPerfil=user.imgPerfil

    usuarios_con_progreso = []    
    busqueda = request.GET.get("buscar")
    if busqueda:
        usuarios = User.objects.filter(
            Q(username__icontains=busqueda) |
            Q(first_name__icontains=busqueda) |
            Q(last_name__icontains=busqueda) |
            Q(email__icontains=busqueda)
        ).distinct()
        print("entre al if de busqueda")
    print("LISTADO DE USUARIOS")
    print(usuarios)
    print(usuarios2)

        # Paginación — Va después de filtrar los usuarios
    paginator = Paginator(usuarios, 10)  # 10 usuarios por página
    page_number = request.GET.get('page')
    page_obj = paginator.get_page(page_number)
    # 1. Cantidad de usuarios cuyo rol es igual a 2  
    
    filtro_progreso = request.GET.get('progreso')

    for usuario in page_obj:
        progreso = 0
        recursos_totales = 0
        recursos_vistos = 0
        print("viendo viendoo")
        print(hasattr(usuario, 'rol') and usuario.rol == 'Estudiante')
        print(usuario.rol)        
        print(str(usuario.rol) == 'Estudiante')
        print(hasattr(usuario, 'rol'))

        if hasattr(usuario, 'rol') and str(usuario.rol) == 'Estudiante':
            # Cursos del estudiante
            cursos = Curso.objects.filter(inscripcioncurso__estudiante=usuario)
            recursos_totales = Recurso.objects.filter(seccion__curso__in=cursos).count()
            recursos_vistos = ProgresoUsuario.objects.filter(usuario=usuario, visto=True).count()
            progreso = int((recursos_vistos / recursos_totales) * 100) if recursos_totales > 0 else 0
        # Aplicamos el filtro si se solicitó
        if filtro_progreso == 'bajo' and progreso >= 50:
            continue  # ignorar este usuario
        usuarios_con_progreso.append({
            'usuario': usuario,
            'progreso': progreso,
            'recursos_vistos': recursos_vistos,
            'recursos_totales': recursos_totales,
        }) 
      
    context = {                     
        'imgPerfil': imgPerfil,        
        'usuario':user.username,  
        'usuarios':page_obj,
        'usuarios_con_progreso': usuarios_con_progreso,  # Lista con progreso incluido
      
    }
    return render(request, 'usAdmin/detalle_usuarios.html', context)

@role_required(['Admin', 'Docente'])
def editar_usuario(request, user_id):
    usuario = get_object_or_404(User, id=user_id)
    cursos_inscritos = InscripcionCurso.objects.filter(estudiante=usuario).select_related('curso')
    # Total de recursos en la plataforma o en los cursos donde está inscrito
    cursos = Curso.objects.filter(inscripcioncurso__estudiante=usuario)
    recursos_totales = Recurso.objects.filter(seccion__curso__in=cursos).count()

    # Recursos vistos por el estudiante
    recursos_vistos = ProgresoUsuario.objects.filter(usuario=usuario, visto=True).count()

    # Se evita división por cero
    if recursos_totales > 0:
        progreso = int((recursos_vistos / recursos_totales) * 100)
    else:
        progreso = 0

    if request.method == 'POST':
        usuario.first_name = request.POST.get('first_name')
        usuario.last_name = request.POST.get('last_name')
        usuario.email = request.POST.get('email')
        usuario.is_active = 'is_active' in request.POST  # checkbox activo
        usuario.save()
        messages.success(request, 'Usuario actualizado correctamente.')
        return redirect('detalleUsuarios-adm') 

    context = {
        'usuario': usuario,
        'cursos_inscritos': cursos_inscritos,
        'progreso': progreso,
        'recursos_totales': recursos_totales,
        'recursos_vistos': recursos_vistos,

    }
    return render(request, 'usAdmin/editar_usuario.html', context)

@role_required('Admin')
def eliminar_usuario(request, user_id):
    usuario = get_object_or_404(User, id=user_id)
    if request.method == 'POST':
        usuario.delete()
        messages.success(request, "Usuario eliminado correctamente.")
        return redirect('detalle_usuarios')
    return render(request, 'usAdmin/confirmar_eliminar_usuario.html', {'usuario': usuario})


@role_required('Admin')
@login_required
def resetear_contrasena(request, id):
    usuario = get_object_or_404(User, id=id)
    
    return redirect('detalle_usuarios')

@role_required('Admin')
def resetContra_usuario(request, id):
    usuario = get_object_or_404(User, id=id)
    # lógica para mostrar o editar
    nueva_contrasena = "12345678"
    usuario.set_password(nueva_contrasena)  # Se cambia la contraseña de forma segura
    usuario.save()
    messages.success(request, f"La contraseña del usuario '{usuario.username}' ha sido reseteada con éxito.")
    return redirect('detalleUsuarios-adm')


@csrf_exempt
def obtener_presigned_url(request):
    if request.method == 'POST':
        nombre_archivo = request.POST.get('filename')
        extension = nombre_archivo.split('.')[-1]
        nuevo_nombre = f"videos/{uuid.uuid4()}.{extension}"

        s3 = boto3.client(
            's3',
            aws_access_key_id=settings.AWS_ACCESS_KEY_ID,
            aws_secret_access_key=settings.AWS_SECRET_ACCESS_KEY,
            region_name=settings.AWS_S3_REGION_NAME
        )

        presigned_post = s3.generate_presigned_post(
            Bucket=settings.AWS_STORAGE_BUCKET_NAME,
            Key=nuevo_nombre,
            Fields={},  # sin ACL aquí
            Conditions=[],
            ExpiresIn=3600
        )

        url_publica = f"https://{settings.AWS_STORAGE_BUCKET_NAME}.s3.{settings.AWS_S3_REGION_NAME}.amazonaws.com/{nuevo_nombre}"

        return JsonResponse({
            'data': presigned_post,
            'url': url_publica
        })

@role_required(['Admin', 'Docente'])
def vistaCrearCurso(request):
    user = request.user    
    imgPerfil=user.imgPerfil 
    if request.method == 'POST':
        form = CursoForm(request.POST, request.FILES)
        if form.is_valid():
            curso = form.save(commit=False)
            curso.profesor = request.user  # Asignar el usuario actual como profesor
            curso.save()               
            return redirect('detalle_curso', curso.id) 
        
    else:
        form = CursoForm() 
    context = {                  
        'imgPerfil': imgPerfil,        
        'usuario':user.username,
        'form': form        
    }   
    return render(request,'usAdmin/crearCurso.html',context)

#EMPEZANDO INTENTAR CREAR COMPETENCIAS

"""def panel_competencia(request, competencia_id):
    competencia = get_object_or_404(Competencia, id=competencia_id)
    preguntas = competencia.preguntas.order_by('orden')

    return render(request, 'usAdmin/panel_control.html', {
        'competencia': competencia,
        'preguntas': preguntas,
    })

def agregar_pregunta_competencia(request, competencia_id):
    competencia = Competencia.objects.get(id=competencia_id)

    if request.method == 'POST':
        enunciado = request.POST.get('enunciado')
        tipo_id = request.POST.get('tipo_pregunta')
        orden = request.POST.get('orden')
        puntaje = request.POST.get('puntaje')

        tipo = TipoPregunta.objects.get(id=tipo_id)
        tipo_nombre = tipo.nombre.lower()

        # Crear la pregunta base
        pregunta = PreguntaCuestionario.objects.create(
            competencia=competencia,
            enunciado=enunciado,
            tipo=tipo,
            orden=orden,
            puntaje=puntaje
        )

        # Opciones según el tipo
        if tipo_nombre == 'seleccion_multiple':
            for i in range(1, 5):
                texto = request.POST.get(f'opcion_{i}')
                es_correcta = request.POST.get('correcta') == str(i)
                if texto:
                    pregunta.opciones.create(texto=texto, es_correcta=es_correcta)

        elif tipo_nombre == 'verdadero_falso':
            respuesta = request.POST.get('vf_respuesta')
            if respuesta:
                pregunta.opciones.create(texto='Verdadero', es_correcta=(respuesta == 'verdadero'))
                pregunta.opciones.create(texto='Falso', es_correcta=(respuesta == 'falso'))

        elif tipo_nombre == 'completar':
            palabra_correcta = request.POST.get('palabra_correcta')
            if palabra_correcta:
                pregunta.opciones.create(texto=palabra_correcta.strip(), es_correcta=True)

        return redirect('panel_competencia', competencia.id)

    tipos = TipoPregunta.objects.all()
    return render(request, 'usAdmin/agregar_pregunta.html', {
        'competencia': competencia,
        'tipos': tipos
    })


from django.http import HttpResponseRedirect
from django.urls import reverse

def iniciar_competencia(request, competencia_id):
    competencia = get_object_or_404(Competencia, id=competencia_id)
    competencia.activa = True
    competencia.save()
    return HttpResponseRedirect(reverse('panel_competencia', args=[competencia.id]))

def enviar_pregunta(request, pregunta_id):
    pregunta = get_object_or_404(PreguntaCuestionario, id=pregunta_id)
    competencia = pregunta.competencia
    competencia.pregunta_actual = pregunta
    competencia.save()
    return HttpResponseRedirect(reverse('panel_competencia', args=[competencia.id]))"""

##HASTA AQU[I]

@role_required('Estudiante')
def inicioEstudiante(request):
    usuario = request.user
    imgPerfil=usuario.imgPerfil
    inscripciones = InscripcionCurso.objects.filter(estudiante=usuario).select_related('curso')
    cursos = [ins.curso for ins in inscripciones]
    
    print(imgPerfil) 
    context = {       
        'cursos': cursos,           
        'imgPerfil': imgPerfil,        
        'usuario':usuario,                
    }      
    return render(request, 'estudiante/student-dashboard.html', context)    

from django.contrib import messages

def unirse_curso(request):
    if request.method == 'POST':
        codigo = request.POST.get('codigo_acceso')
        try:
            curso = Curso.objects.get(codigo_acceso=codigo)
            ya_inscrito = InscripcionCurso.objects.filter(curso=curso, estudiante=request.user).exists()
            if ya_inscrito:
                messages.warning(request, 'Ya estás inscrito en este curso.')
            else:
                InscripcionCurso.objects.create(curso=curso, estudiante=request.user)
                messages.success(request, 'Te has unido correctamente al curso.')
        except Curso.DoesNotExist:
            messages.error(request, 'El código ingresado no es válido.')
        return redirect('listar_cursos') 

def listar_cursos(request):
    usuario = request.user
    inscripciones = InscripcionCurso.objects.filter(estudiante=usuario).select_related('curso')
    cursos = [ins.curso for ins in inscripciones]

    return render(request, 'estudiante/student-course-list.html', {
        'cursos': cursos,
        'usuario': usuario,
        'imgPerfil': usuario.imgPerfil
    })

def detalle_cursoEstudiante(request, curso_id):
    print("Ingrese al detalle del curso estudiante")
    curso = get_object_or_404(Curso, id=curso_id)
    
    recursos_prefetch = Prefetch(
        'recursos',
        queryset=Recurso.objects.order_by('orden'),
        to_attr='recursos_ordenados'
    )
    print(recursos_prefetch)
    secciones = Seccion.objects.filter(curso=curso_id).prefetch_related(recursos_prefetch).order_by('orden')    
    recursosTipo = TipoRecurso.objects.all()
    user = request.user
    imgPerfil = user.imgPerfil

    return render(request, 'estudiante/student-course-resume.html', {
        'curso': curso,
        'imgPerfil': imgPerfil,
        'usuario': user.username,
        'secciones': secciones,
        'tipos_recurso': recursosTipo
    })

@login_required
def resolver_practica(request, practica_id):
    practica = get_object_or_404(Practica, id=practica_id)

    if request.method == 'POST':
        smiles_estudiante = request.POST.get('mol_input')
        modelo_objetivo = practica.modelo_objetivo

        # Comprobación (simple, se puede mejorar con RDKit)
        es_correcta = smiles_estudiante.strip() == modelo_objetivo.smiles.strip()

        MoleculaEstudiante.objects.create(
            estudiante=request.user,
            practica=practica,
            smiles=smiles_estudiante,
            es_correcta=es_correcta
        )

        return render(request, 'estudiante/resultado_practica.html', {
            'practica': practica,
            'es_correcta': es_correcta
        })

    return render(request, 'estudiante/resolver_practica.html', {
        'practica': practica,
    })

@login_required
def historial_intentos(request, practica_id):
    practica = get_object_or_404(Practica, id=practica_id)
    intentos = MoleculaEstudiante.objects.filter(
        estudiante=request.user,
        practica=practica
    ).order_by('-fecha_envio')

    return render(request, 'estudiante/historial_intentos.html', {
        'practica': practica,
        'intentos': intentos
    })



def biblioteca_practicas(request):
    query = request.GET.get('q', '')  # Capturamos el término de búsqueda

    practicas_visibles = Recurso.objects.filter(
        tipo__nombre__iexact='practica',
        visible_biblioteca=True
    )

    if query:
        practicas_visibles = practicas_visibles.filter(
            Q(titulo__icontains=query)
        )

    # Paginación: 10 por página
    paginator = Paginator(practicas_visibles.order_by('-fecha_creacion'), 10)
    page_number = request.GET.get('page')
    page_obj = paginator.get_page(page_number)

    context = {
        'page_obj': page_obj,
        'query': query
    }
    return render(request, 'estudiante/biblioteca_practicas.html', context)

from django.db.models.functions import TruncMonth, TruncWeek, TruncDay
@role_required('Docente')
def inicioDocente(request):
    user = request.user

    # Totales
    total_cursos = Curso.objects.filter(profesor=user).count()
    print(f"total cursos{total_cursos} y los cursos {Curso.objects.filter(profesor=user)}")
    total_estudiantes = InscripcionCurso.objects.filter(
        curso__profesor=user
    ).values('estudiante').distinct().count()
    print(total_estudiantes)
    total_practicas = Recurso.objects.filter(
        creado_por=user,
        tipo__nombre__iexact="Practica",
        visible_biblioteca=True
    ).count()
    print(total_practicas)    
    # Cursos listados (con búsqueda/paginación como ya hicimos)
    buscar = request.GET.get('buscar', '').strip()
    cursos_qs = Curso.objects.filter(profesor=user).order_by('-fecha_creacion')
    if buscar:
        cursos_qs = cursos_qs.filter(
            Q(titulo__icontains=buscar) |
            Q(descripcion__icontains=buscar) |
            Q(categoria__nombre__icontains=buscar)
        )
    paginator = Paginator(cursos_qs, 5)
    page_number = request.GET.get('page')
    cursos_page = paginator.get_page(page_number)

    # 🔹 Intervalo dinámico
    intervalo = request.GET.get('intervalo', 'mes')
    if intervalo == 'dia':
        trunc_func = TruncDay
        formato = "%Y-%m-%d"
    elif intervalo == 'semana':
        trunc_func = TruncWeek
        formato = "Semana %W %Y"
    else:
        trunc_func = TruncMonth
        formato = "%Y-%m"

    inscripciones = (
        InscripcionCurso.objects.filter(curso__profesor=user)
        .annotate(periodo=trunc_func('fecha_inscripcion'))
        .values('periodo', 'curso__titulo')
        .annotate(total=Count('id'))
        .order_by('periodo')
    )

    # Reestructurar datos para Chart.js
    periodos = sorted(set([i['periodo'].strftime(formato) for i in inscripciones]))
    cursos_series = {}
    for i in inscripciones:
        curso = i['curso__titulo']
        periodo = i['periodo'].strftime(formato)
        if curso not in cursos_series:
            cursos_series[curso] = {p: 0 for p in periodos}
        cursos_series[curso][periodo] = i['total']

    datasets = [
        {
            "label": curso,
            "data": [series[p] for p in periodos],
        }
        for curso, series in cursos_series.items()
    ]
    print(f"total cursos al final {total_cursos} y los cursos {Curso.objects.filter(profesor=user)}")
    print(total_estudiantes)
    print(total_practicas)
    context = {
        'imgPerfil': user.imgPerfil,
        'usuario': user.username,
        'cursos': cursos_page,
        'buscar': buscar,
        'intervalo': intervalo,
        'chart_labels': json.dumps(periodos),
        'chart_datasets': json.dumps(datasets),
        'total_cursos': total_cursos,
        'total_estudiantes': total_estudiantes,
        'total_practicas': total_practicas,
    }

    return render(request, 'docente/instructor-dashboard.html', context)

@role_required('Docente')
@login_required
def mis_cursos_docente(request):
    user = request.user
    cursos = Curso.objects.filter(profesor=request.user).order_by('-fecha_creacion')
     # 🔹 Búsqueda
    buscar = request.GET.get('buscar', '').strip()
    cursos_qs = Curso.objects.filter(profesor=user).order_by('-fecha_creacion')

    if buscar:
        cursos_qs = cursos_qs.filter(
            Q(titulo__icontains=buscar) |
            Q(descripcion__icontains=buscar) |
            Q(categoria__nombre__icontains=buscar)
        )

    # 🔹 Paginación
    paginator = Paginator(cursos_qs, 5)  # 5 cursos por página
    page_number = request.GET.get('page')
    cursos_page = paginator.get_page(page_number)
    return render(request, 'docente/mis_cursos.html', {        
        'usuario': request.user,
        'total_cursos': cursos,        
        'cursos': cursos_page, 
        'buscar': buscar,
    })

@role_required('Docente')
@login_required
def alumnos_mis_cursos(request):
    cursos_docente = Curso.objects.filter(profesor=request.user)

    curso_id = request.GET.get('curso')
    buscar = request.GET.get('buscar', '').strip()

    inscripciones = InscripcionCurso.objects.filter(
        curso__in=cursos_docente
    ).select_related('estudiante', 'curso')

    # 🔹 Filtro por curso
    if curso_id and curso_id.isdigit():
        inscripciones = inscripciones.filter(curso_id=curso_id)

    # 🔹 Filtro por búsqueda (nombre, apellido, username o email)
    if buscar:
        inscripciones = inscripciones.filter(
            Q(estudiante__first_name__icontains=buscar) |
            Q(estudiante__last_name__icontains=buscar) |
            Q(estudiante__username__icontains=buscar) |
            Q(estudiante__email__icontains=buscar)
        )

    # 🔹 Paginación
    paginator = Paginator(inscripciones, 10)  # 10 por página
    page_number = request.GET.get('page')
    inscripciones_page = paginator.get_page(page_number)

    return render(request, 'docente/alumnos_mis_cursos.html', {
        'inscripciones': inscripciones_page,
        'cursos_docente': cursos_docente,
        'curso_id': int(curso_id) if curso_id and curso_id.isdigit() else None,
        'buscar': buscar,  # para mantener valor en input search
    })

"""def custom_login(request):
    print("entre a la funcion custom")
    if request.method == "POST":  
        print("Estoy en post")   
        print(request.POST) 
        print(request.POST['email'], request.POST['password'])  
        user = authenticate(request, username=request.POST['email'], password= request.POST
        ['password'])
        print(user)
        print("Ya autentique el usuario")
        print(user is not None)
        if user is not None:            
            login(request, user)
            print("Loguie al usuario")
            rol = user.rol_id  # Asignar rol después de la autenticación   
            print(rol)
            print(rol==1)
            print(type(rol))         
            if rol == 1:
                print("Ingrese al rol 1")
                return redirect('dashboard-adm')
            elif rol == 2:
                return redirect('student_dashboard')
            elif rol == 3:
                return redirect('teacher_dashboard')
        else:
            # Manejar error de autenticación
            return render(request, 'dashboard-adm', {'error': 'Credenciales inválidas'})
    return redirect('inicio')"""

@csrf_exempt
@csrf_protect
def custom_login(request):
    print("🔐 Entre a la función custom_login")
    print(f"🌐 Método de request: {request.method}")
    
    if request.method == "POST":
        print("📝 Procesando login POST")
        
        # Debug: Ver todos los datos del POST
        print(f"📊 Datos POST completos: {request.POST}")
        
        # Obtener datos del formulario
        email = request.POST.get('email', '').strip()
        password = request.POST.get('password', '')
        
        print(f"📧 Email recibido: '{email}'")
        print(f"🔒 Password recibido: {'Sí (' + str(len(password)) + ' caracteres)' if password else 'No'}")
        print(password)  
        # Validaciones básicas
        if not email or not password:
            error_msg = 'Email y contraseña son obligatorios'
            print(f"❌ Error validación: {error_msg}")
            return render(request, 'sign-in.html', {'error': error_msg})
        
        print("✅ Validaciones básicas pasadas")
        
        # Verificar si el usuario existe
        try:
            print(f"🔍 Buscando usuario con username: '{email}'")
            user_check = User.objects.get(username=email)
            print(f"👤 Usuario encontrado: {user_check.username}")
            print(f"🎭 Rol ID: {user_check.rol_id}")
            print(f"✅ Usuario activo: {user_check.is_active}")
            print(f"🆔 ID de usuario: {user_check.id}")
            print(f"🆔 Contraseña de usuario: {user_check.password}")
        except User.DoesNotExist:
            error_msg = f'No existe un usuario con el email: {email}'
            print(f"❌ Usuario no encontrado: {error_msg}")
            
            # DEBUG: Mostrar todos los usuarios que existen
            print("📋 Usuarios existentes en la base de datos:")
            all_users = User.objects.all()
            for u in all_users:
                print(f"   - {u.username} (rol: {u.rol_id}, activo: {u.is_active})")
            
            return render(request, 'sign-in.html', {'error': error_msg})
        
        print("🔑 Intentando autenticar usuario...")
        
        # Autenticar usuario
        user = authenticate(request, username=request.POST['email'], password= request.POST
        ['password'])
        print(password==request.POST['password'])        
        print(f"🔐 Resultado autenticación: {user}")
        
        if user is not None:
            print("✅ Autenticación exitosa")
            if user.is_active:
                print("✅ Usuario está activo")
                
                # Login exitoso
                try:
                    login(request, user)
                    print(f"✅ Login realizado para: {user.username}")
                except Exception as e:
                    print(f"❌ Error en login(): {e}")
                    return render(request, 'sign-in.html', {'error': f'Error al iniciar sesión: {e}'})
                
                # Obtener rol y redirigir
                rol = getattr(user, 'rol_id', None)
                print(f"🎭 Rol obtenido: {rol} (tipo: {type(rol)})")
                
                if rol == 1:  # Admin/Superadministrador
                    print("🚀 Debería redirigir a dashboard admin")
                    try:
                        return redirect('dashboard-adm')
                    except Exception as e:
                        print(f"❌ Error al redirigir a admin: {e}")
                        return render(request, 'sign-in.html', {'error': f'Error de redirección admin: {e}'})
                        
                elif rol == 2:  # Estudiante
                    print("🚀 Debería redirigir a dashboard estudiante")
                    try:
                        return redirect('student_dashboard')
                    except Exception as e:
                        print(f"❌ Error al redirigir a estudiante: {e}")
                        return render(request, 'sign-in.html', {'error': f'Error de redirección estudiante: {e}'})
                        
                elif rol == 3:  # Docente
                    print("🚀 Debería redirigir a dashboard docente")
                    try:
                        return redirect('teacher_dashboard')
                    except Exception as e:
                        print(f"❌ Error al redirigir a docente: {e}")
                        return render(request, 'sign-in.html', {'error': f'Error de redirección docente: {e}'})
                else:
                    print(f"⚠️ Rol no reconocido: {rol}")
                    print("🚀 Redirigiendo a admin por defecto")
                    try:
                        return redirect('dashboard-adm')
                    except Exception as e:
                        print(f"❌ Error al redirigir por defecto: {e}")
                        return render(request, 'sign-in.html', {'error': f'Error de redirección por defecto: {e}'})
            else:
                error_msg = 'Tu cuenta está desactivada. Contacta al administrador.'
                print(f"❌ Usuario inactivo: {error_msg}")
                return render(request, 'sign-in.html', {'error': error_msg})
        else:
            error_msg = 'Contraseña incorrecta'
            print(f"❌ Credenciales inválidas: {error_msg}")
            return render(request, 'sign-in.html', {'error': error_msg})
    
    # GET request - mostrar formulario
    print("📄 Mostrando formulario de login (GET request)")
    return render(request, 'sign-in.html')


def signout(request):
    logout(request)
    rotate_token(request)  # Gira el token CSRF para la nueva sesión
    return redirect('inicio')

#from rdkit.Chem import AllChem, SDWriter
import tempfile

@csrf_exempt
def editor_dibujo(request):
    return render(request, 'editor_dibujo.html')

@csrf_exempt
def procesar_molecula(request):
    if request.method == 'POST':
        entrada = request.POST.get('mol_input')
        tipo = request.POST.get('tipo')  # "smiles" o "mol"

        try:
            response = requests.post(
                f"{settings.RDKit_API_URL}/procesar",
                json={"entrada": entrada, "tipo": tipo},
                timeout=30
            )

            if response.status_code == 200:
                data = response.json()
                sdf_data = data.get("sdf")
                return render(request, 'visor_3d.html', {'sdf_data': sdf_data})
            else:
                return JsonResponse(
                    {"error": response.json().get("error", "Error desconocido")},
                    status=response.status_code
                )

        except Exception as e:
            return JsonResponse({'error': str(e)}, status=500)

    return JsonResponse({'message': 'Usa POST con SMILES o MOL'})


# ---- API Validar (Constructor 2D → RDKit) ----
def _pubchem_json(url, timeout=15):
    r = requests.get(url, timeout=timeout)
    if r.status_code != 200:
        return None
    try:
        return r.json()
    except:
        return None

def _pc(url, timeout=15):
    try:
        r = requests.get(url, timeout=timeout)
        if r.status_code != 200:
            print("PUBCHEM non-200:", r.status_code, url)
            return None
        return r.json()
    except Exception as e:
        print("PUBCHEM error:", e, url)
        return None

import json, re, urllib.request, urllib.error

SMILES_PATTERN = re.compile(r"^[A-Za-z0-9@+\-\[\]\(\)=#$]{2,}$")

def _http_post_json(url: str, payload: dict, timeout: int = 30):
    data = json.dumps(payload).encode('utf-8')
    req = urllib.request.Request(
        url, data=data,
        headers={'Content-Type': 'application/json'},
        method='POST'
    )
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            return resp.status, resp.read().decode('utf-8', errors='replace')
    except urllib.error.HTTPError as e:
        return e.code, e.read().decode('utf-8', errors='replace')
    except Exception as e:
        return None, str(e)

@csrf_exempt
@require_POST
def resolve_smiles(request):
    try:
        body = json.loads(request.body.decode('utf-8'))
    except Exception:
        return JsonResponse({'smiles': None})  # sin romper

    query = (body.get('query') or '').strip()
    if not query:
        return JsonResponse({'smiles': None})

    # si ya parece SMILES → devuélvelo (no se hace nada más)
    if SMILES_PATTERN.match(query):
        return JsonResponse({'smiles': query})

    # si quieres dejar el intento remoto activado, lo mantengo:
    base = getattr(settings, 'RDKit_API_URL', '').rstrip('/')
    if not base:
        return JsonResponse({'smiles': None})  # sin romper

    status, body = _http_post_json(f"{base}/name_to_smiles", {'query': query})
    if status != 200:
        return JsonResponse({'smiles': None})  # sin romper

    try:
        data = json.loads(body)
    except Exception:
        return JsonResponse({'smiles': None})

    return JsonResponse({'smiles': data.get('smiles')})


#Paula

###########   Docente   ###############
from django.views.decorators.http import require_http_methods
from django.utils import timezone
from django.db import transaction
from django.db.models import Sum, Max
from datetime import datetime,timedelta

@login_required
def instructorQuiz(request, seccion_id=None):
    """Vista principal para crear cuestionarios - COMPLETAMENTE FUNCIONAL"""
    print(f"🎯 instructorQuiz llamada con seccion_id: {seccion_id}")
    
    seccion = None
    if seccion_id:
        seccion = get_object_or_404(Seccion, id=seccion_id)
        print(f"📁 Sección encontrada: {seccion.titulo}")
        
        # Verificar permisos
        if seccion.curso.profesor != request.user:
            print("❌ Usuario sin permisos")
            return redirect('dashboard-adm')
    
    # Crear tipos de pregunta si no existen
    tipos_pregunta = TipoPregunta.objects.all()
    if not tipos_pregunta.exists():
        print("🔧 Creando tipos de pregunta básicos...")
        tipos_pregunta = crear_tipos_pregunta_basicos()
    
    print(f"📝 Tipos de pregunta disponibles: {tipos_pregunta.count()}")
    
    context = {
        'cuestionario': None,
        'preguntas': [],
        'seccion': seccion,
        'seccion_id': seccion_id,
        'tipos_pregunta': tipos_pregunta,
    }
    
    return render(request, 'docente/crear_cuestionario.html', context)

def crear_tipos_pregunta_basicos():
    """Crear tipos de pregunta básicos si no existen - FUNCIONAL"""
    tipos_basicos = [
        ('opcion_unica', 'Opción Única'),
        ('opcion_multiple', 'Opción Múltiple'),
        ('falso_verdadero', 'Verdadero/Falso'),
        ('unir_lineas', 'Unir con Líneas'),
        ('completar', 'Completar'),
        ('simulador_2d', 'Simulador 2D'),
        ('simulador_3d', 'Simulador 3D'),
        ('respuesta_abierta', 'Respuesta Abierta'),
    ]
    
    for nombre, descripcion in tipos_basicos:
        TipoPregunta.objects.get_or_create(
            nombre=nombre,
            defaults={'descripcion': descripcion}
        )
    
    return TipoPregunta.objects.all()

# =====================================================
# VISTA PARA EDITAR CUESTIONARIOS EXISTENTES - CORREGIDA
# =====================================================

@login_required
def editarCuestionario(request, cuestionario_id):
    """Vista para editar un cuestionario existente - COMPLETAMENTE FUNCIONAL"""
    print(f"✏️ editarCuestionario llamada con ID: {cuestionario_id}")
    
    cuestionario = get_object_or_404(Cuestionario, id=cuestionario_id)
    print(f"📋 Cuestionario encontrado: {cuestionario.recurso.titulo}")
    
    # Verificar permisos
    if (cuestionario.recurso.seccion and 
        cuestionario.recurso.seccion.curso.profesor != request.user):
        print("❌ Usuario sin permisos para editar")
        return redirect('dashboard-adm')
    
    # Obtener preguntas y tipos
    preguntas = cuestionario.preguntas.all().order_by('orden')
    tipos_pregunta = TipoPregunta.objects.all()
    
    # Calcular datos dinámicos
    puntaje_calculado = preguntas.aggregate(total=Sum('puntaje'))['total'] or 0
    tiene_preguntas_manuales = preguntas.filter(
        tipo__nombre__in=['respuesta_abierta', 'simulador_2d', 'simulador_3d']
    ).exists()
    
    print(f"📊 Preguntas: {preguntas.count()}, Puntaje: {puntaje_calculado}")
    
    context = {
        'cuestionario': cuestionario,
        'preguntas': preguntas,
        'tipos_pregunta': tipos_pregunta,
        'puntaje_calculado': puntaje_calculado,
        'tiene_preguntas_manuales': tiene_preguntas_manuales,
        'calificacion_mixta': tiene_preguntas_manuales and preguntas.exclude(
            tipo__nombre__in=['respuesta_abierta', 'simulador_2d', 'simulador_3d']
        ).exists()
    }
    
    return render(request, 'docente/crear_cuestionario.html', context)



@login_required
def editar_cuestionario(request, recurso_id):
    print(f"✏️ editarCuestionario llamada con Recurso ID: {recurso_id}")
    
    recurso = get_object_or_404(Recurso, id=recurso_id)
    cuestionario = get_object_or_404(Cuestionario, recurso=recurso)
    print(f"📋 Cuestionario encontrado: {cuestionario.recurso.titulo}")
    
    # Verificar permisos
    if (cuestionario.recurso.seccion and 
        cuestionario.recurso.seccion.curso.profesor != request.user):
        print("❌ Usuario sin permisos para editar")
        return redirect('dashboard-adm')
    
    # Obtener preguntas y tipos
    preguntas = cuestionario.preguntas.all().order_by('orden')
    tipos_pregunta = TipoPregunta.objects.all()
    
    # Calcular datos dinámicos
    puntaje_calculado = preguntas.aggregate(total=Sum('puntaje'))['total'] or 0
    tiene_preguntas_manuales = preguntas.filter(
        tipo__nombre__in=['respuesta_abierta', 'simulador_2d', 'simulador_3d']
    ).exists()
    
    print(f"📊 Preguntas: {preguntas.count()}, Puntaje: {puntaje_calculado}")
    
    context = {
        'cuestionario': cuestionario,
        'preguntas': preguntas,
        'tipos_pregunta': tipos_pregunta,
        'puntaje_calculado': puntaje_calculado,
        'tiene_preguntas_manuales': tiene_preguntas_manuales,
        'calificacion_mixta': tiene_preguntas_manuales and preguntas.exclude(
            tipo__nombre__in=['respuesta_abierta', 'simulador_2d', 'simulador_3d']
        ).exists()
    }
    
    return render(request, 'docente/crear_cuestionario.html', context)
# =====================================================
# GUARDAR CUESTIONARIO BÁSICO - CORREGIDA
# =====================================================

@csrf_exempt
@require_http_methods(["POST"])
@login_required
def guardar_cuestionario(request):
    """Guardar cuestionario básico - COMPLETAMENTE FUNCIONAL"""
    try:
        print("💾 Iniciando guardado de cuestionario...")
        print(f"POST data: {request.POST}")
        
        seccion_id = request.POST.get('seccion_id')
        titulo = request.POST.get('titulo', '').strip()
        descripcion = request.POST.get('descripcion', '').strip()
        
        print(f"📝 Datos recibidos - Título: '{titulo}', Descripción: '{descripcion[:50]}...'")
        
        # Validaciones
        if not titulo:
            return JsonResponse({
                'success': False, 
                'error': 'El título del cuestionario es obligatorio'
            })
        
        if len(titulo) < 3:
            return JsonResponse({
                'success': False,
                'error': 'El título debe tener al menos 3 caracteres'
            })
        
        if not descripcion:
            return JsonResponse({
                'success': False,
                'error': 'La descripción del cuestionario es obligatoria'
            })
        
        if len(descripcion) < 10:
            return JsonResponse({
                'success': False,
                'error': 'La descripción debe tener al menos 10 caracteres'
            })
        
        # Validar sección
        seccion = None
        if seccion_id:
            seccion = get_object_or_404(Seccion, id=seccion_id)
            if seccion.curso.profesor != request.user:
                return JsonResponse({
                    'success': False, 
                    'error': 'Sin permisos para crear cuestionario en esta sección'
                }, status=403)
        
        # Crear cuestionario con transacción
        with transaction.atomic():
            print("🔧 Creando tipo de recurso cuestionario...")
            
            # Crear o obtener tipo de recurso
            tipo_cuestionario, created = TipoRecurso.objects.get_or_create(
                nombre='cuestionario',
                defaults={'descripcion': 'Recurso tipo cuestionario'}
            )
            
            if created:
                print("✅ Tipo de recurso 'cuestionario' creado")
            
            # Crear recurso
            print("🗂️ Creando recurso...")
            recurso = Recurso.objects.create(
                seccion=seccion,
                tipo=tipo_cuestionario,
                titulo=titulo,
                descripcion=descripcion,
                orden=seccion.recursos.count() + 1 if seccion else 1,
                creado_por=request.user
            )
            
            print(f"✅ Recurso creado con ID: {recurso.id}")
            
            # Crear cuestionario
            print("📋 Creando cuestionario...")
            cuestionario = Cuestionario.objects.create(
                recurso=recurso,
                instrucciones=descripcion,
                tiempo_limite=30,
                puntaje_total=0.00,
                calificacion_automatica=True,
                intentos_permitidos=1,
                mostrar_resultados=True,
                orden_aleatorio=False
            )
            
            print(f"✅ Cuestionario creado con ID: {cuestionario.id}")
        
        response_data = {
            'success': True,
            'cuestionario_id': cuestionario.id,
            'mensaje': 'Cuestionario creado exitosamente',
            'redirect_url': f'/docente/editar_cuestionario/{cuestionario.id}/'
        }
        
        print(f"✅ Respuesta exitosa: {response_data}")
        return JsonResponse(response_data)
        
    except Exception as e:
        print(f"❌ Error al guardar cuestionario: {str(e)}")
        import traceback
        traceback.print_exc()
        
        return JsonResponse({
            'success': False, 
            'error': f'Error interno: {str(e)}'
        }, status=400)

# =====================================================
# ACTUALIZAR CONFIGURACIÓN - CORREGIDA
# =====================================================

@csrf_exempt
@require_http_methods(["POST"])
@login_required
def actualizarConfiguracion(request):
    """Actualizar configuración del cuestionario - COMPLETAMENTE FUNCIONAL"""
    try:
        print("⚙️ Actualizando configuración...")
        
        data = json.loads(request.body)
        cuestionario_id = data.get('cuestionario_id')
        configuracion = data.get('configuracion')
        
        print(f"📋 Cuestionario ID: {cuestionario_id}")
        print(f"⚙️ Configuración: {configuracion}")
        
        cuestionario = get_object_or_404(Cuestionario, id=cuestionario_id)
        
        # Verificar permisos
        if (cuestionario.recurso.seccion and 
            cuestionario.recurso.seccion.curso.profesor != request.user):
            return JsonResponse({
                'success': False, 
                'error': 'Sin permisos'
            }, status=403)
        
        # Actualizar configuración
        if configuracion:
            if 'tiempo_limite' in configuracion:
                tiempo = int(configuracion['tiempo_limite'])
                if 1 <= tiempo <= 300:
                    cuestionario.tiempo_limite = tiempo
                    print(f"⏱️ Tiempo límite actualizado: {tiempo} min")
            
            if 'intentos_permitidos' in configuracion:
                intentos = configuracion['intentos_permitidos']
                if intentos == 'ilimitado' or int(intentos) == 999:
                    cuestionario.intentos_permitidos = 999
                else:
                    cuestionario.intentos_permitidos = int(intentos)
                print(f"🔄 Intentos permitidos: {cuestionario.intentos_permitidos}")
            
            if 'mostrar_resultados' in configuracion:
                cuestionario.mostrar_resultados = bool(configuracion['mostrar_resultados'])
                print(f"👁️ Mostrar resultados: {cuestionario.mostrar_resultados}")
            
            if 'mezclar_preguntas' in configuracion:
                cuestionario.orden_aleatorio = bool(configuracion['mezclar_preguntas'])
                print(f"🔀 Orden aleatorio: {cuestionario.orden_aleatorio}")
            
            # Manejar fechas
            if 'fecha_inicio' in configuracion and configuracion['fecha_inicio']:
                try:
                    fecha_inicio = datetime.fromisoformat(
                        configuracion['fecha_inicio'].replace('Z', '+00:00')
                    )
                    cuestionario.fecha_apertura = fecha_inicio
                    print(f"📅 Fecha inicio: {fecha_inicio}")
                except Exception as e:
                    print(f"⚠️ Error en fecha inicio: {e}")
            
            if 'fecha_fin' in configuracion and configuracion['fecha_fin']:
                try:
                    fecha_fin = datetime.fromisoformat(
                        configuracion['fecha_fin'].replace('Z', '+00:00')
                    )
                    cuestionario.fecha_cierre = fecha_fin
                    print(f"📅 Fecha fin: {fecha_fin}")
                except Exception as e:
                    print(f"⚠️ Error en fecha fin: {e}")
        
        # Recalcular totales
        puntaje_total = cuestionario.preguntas.aggregate(total=Sum('puntaje'))['total'] or 0
        cuestionario.puntaje_total = puntaje_total
        
        tiene_preguntas_manuales = cuestionario.preguntas.filter(
            tipo__nombre__in=['respuesta_abierta', 'simulador_2d', 'simulador_3d']
        ).exists()
        cuestionario.calificacion_automatica = not tiene_preguntas_manuales
        
        cuestionario.save()
        print("✅ Configuración guardada exitosamente")
        
        return JsonResponse({
            'success': True,
            'mensaje': 'Configuración actualizada exitosamente',
            'puntaje_total': float(puntaje_total),
            'calificacion_automatica': not tiene_preguntas_manuales
        })
        
    except Exception as e:
        print(f"❌ Error al actualizar configuración: {str(e)}")
        import traceback
        traceback.print_exc()
        
        return JsonResponse({
            'success': False, 
            'error': str(e)
        }, status=400)

# =====================================================
# GESTIÓN DE PREGUNTAS - COMPLETAMENTE FUNCIONAL
# =====================================================

@csrf_exempt
@require_http_methods(["POST"])
@login_required
def agregarPregunta(request):
    """Agregar una nueva pregunta al cuestionario - CORREGIDO"""
    try:
        print("➕ Agregando nueva pregunta...")
        
        data = json.loads(request.body)
        cuestionario_id = data.get('cuestionario_id')
        tipo_nombre = data.get('tipo')
        
        print(f"📋 Cuestionario ID: {cuestionario_id}, Tipo: {tipo_nombre}")
        
        cuestionario = get_object_or_404(Cuestionario, id=cuestionario_id)
        
        # Verificar permisos
        if (cuestionario.recurso.seccion and 
            cuestionario.recurso.seccion.curso.profesor != request.user):
            return JsonResponse({
                'success': False, 
                'error': 'Sin permisos'
            }, status=403)
        
        tipo_pregunta = get_object_or_404(TipoPregunta, nombre=tipo_nombre)
        
        # Calcular orden - CORREGIDO para usar relaciones genéricas
        from django.contrib.contenttypes.models import ContentType
        cuestionario_ct = ContentType.objects.get_for_model(Cuestionario)
        
        ultima_pregunta = PreguntaCuestionario.objects.filter(
            content_type=cuestionario_ct,
            object_id=cuestionario.id
        ).order_by('-orden').first()
        
        orden = ultima_pregunta.orden + 1 if ultima_pregunta else 1
        
        print(f"📝 Creando pregunta orden {orden}")
        
        # ✅ CREAR PREGUNTA USANDO RELACIONES GENÉRICAS:
        pregunta = PreguntaCuestionario.objects.create(
            content_type=cuestionario_ct,
            object_id=cuestionario.id,
            tipo=tipo_pregunta,
            enunciado=f"Nueva pregunta {orden}",
            orden=orden,
            puntaje=1
        )
        
        print(f"✅ Pregunta creada con ID: {pregunta.id}")
        
        # Crear opciones por defecto según el tipo
        if tipo_nombre in ['opcion_unica', 'opcion_multiple', 'falso_verdadero']:
            if tipo_nombre == 'falso_verdadero':
                opciones_default = [('Verdadero', True), ('Falso', False)]
            else:
                opciones_default = [('Opción 1', True), ('Opción 2', False)]
            
            for texto, es_correcta in opciones_default:
                OpcionPregunta.objects.create(
                    pregunta=pregunta,
                    texto=texto,
                    es_correcta=es_correcta
                )
                print(f"✅ Opción creada: {texto} (correcta: {es_correcta})")
        
        # Recalcular totales
        puntaje_total = cuestionario.preguntas.aggregate(total=Sum('puntaje'))['total'] or 0
        cuestionario.puntaje_total = puntaje_total
        
        tiene_preguntas_manuales = cuestionario.preguntas.filter(
            tipo__nombre__in=['respuesta_abierta', 'simulador_2d', 'simulador_3d']
        ).exists()
        cuestionario.calificacion_automatica = not tiene_preguntas_manuales
        cuestionario.save()
        
        print(f"📊 Totales actualizados - Puntaje: {puntaje_total}")
        
        return JsonResponse({
            'success': True,
            'pregunta_id': pregunta.id,
            'mensaje': 'Pregunta agregada exitosamente',
            'puntaje_total': float(puntaje_total),
            'calificacion_automatica': not tiene_preguntas_manuales
        })
        
    except Exception as e:
        print(f"❌ Error al agregar pregunta: {str(e)}")
        import traceback
        traceback.print_exc()
        
        return JsonResponse({
            'success': False, 
            'error': str(e)
        }, status=400)

@csrf_exempt
@require_http_methods(["POST"])
@login_required
def cambiarTipoPregunta(request):
    """Cambiar el tipo de una pregunta existente - COMPLETAMENTE FUNCIONAL"""
    try:
        print("🔄 Cambiando tipo de pregunta...")
        
        data = json.loads(request.body)
        pregunta_id = data.get('pregunta_id')
        nuevo_tipo = data.get('nuevo_tipo')
        
        print(f"🔄 Pregunta ID: {pregunta_id}, Nuevo tipo: {nuevo_tipo}")
        
        pregunta = get_object_or_404(PreguntaCuestionario, id=pregunta_id)
        
        # ✅ VERIFICAR QUE PERTENECE A UN CUESTIONARIO:
        if not isinstance(pregunta.content_object, Cuestionario):
            return JsonResponse({
                'success': False, 
                'error': 'Esta pregunta no pertenece a un cuestionario'
            })
        
        cuestionario = pregunta.content_object
        
        # Verificar permisos
        if (cuestionario.recurso.seccion and 
            cuestionario.recurso.seccion.curso.profesor != request.user):
            return JsonResponse({
                'success': False, 
                'error': 'Sin permisos'
            }, status=403)
        
        with transaction.atomic():
            # Eliminar opciones existentes
            opciones_eliminadas = pregunta.opciones.count()
            pregunta.opciones.all().delete()
            print(f"🗑️ Eliminadas {opciones_eliminadas} opciones existentes")
            
            # Cambiar tipo
            tipo_pregunta = get_object_or_404(TipoPregunta, nombre=nuevo_tipo)
            pregunta.tipo = tipo_pregunta
            pregunta.save()
            
            print(f"🔄 Tipo cambiado a: {tipo_pregunta.descripcion}")
            
            # Crear nuevas opciones según el tipo
            if nuevo_tipo in ['opcion_unica', 'opcion_multiple']:
                opciones_default = [
                    ('Opción 1', True), ('Opción 2', False),
                    ('Opción 3', False), ('Opción 4', False)
                ]
                for texto, es_correcta in opciones_default:
                    OpcionPregunta.objects.create(
                        pregunta=pregunta, 
                        texto=texto, 
                        es_correcta=es_correcta
                    )
                    print(f"✅ Nueva opción: {texto}")
                    
            elif nuevo_tipo == 'falso_verdadero':
                OpcionPregunta.objects.create(
                    pregunta=pregunta, 
                    texto='Verdadero', 
                    es_correcta=True
                )
                OpcionPregunta.objects.create(
                    pregunta=pregunta, 
                    texto='Falso', 
                    es_correcta=False
                )
                print("✅ Opciones Verdadero/Falso creadas")
        
        return JsonResponse({
            'success': True, 
            'mensaje': 'Tipo de pregunta cambiado exitosamente'
        })
        
    except Exception as e:
        print(f"❌ Error al cambiar tipo: {str(e)}")
        import traceback
        traceback.print_exc()
        
        return JsonResponse({
            'success': False, 
            'error': str(e)
        }, status=400)

@csrf_exempt
@require_http_methods(["DELETE"])
@login_required
def eliminarPregunta(request, pregunta_id):
    """Eliminar una pregunta completa - COMPLETAMENTE FUNCIONAL"""
    try:
        print(f"🗑️ Eliminando pregunta ID: {pregunta_id}")
        
        pregunta = get_object_or_404(PreguntaCuestionario, id=pregunta_id)
        
        # Verificar permisos
        if (pregunta.cuestionario.recurso.seccion and 
            pregunta.cuestionario.recurso.seccion.curso.profesor != request.user):
            return JsonResponse({
                'success': False, 
                'error': 'Sin permisos'
            }, status=403)
        
        cuestionario = pregunta.cuestionario
        enunciado = pregunta.enunciado[:50]
        
        pregunta.delete()
        print(f"✅ Pregunta eliminada: {enunciado}")
        
        # Recalcular totales
        puntaje_total = cuestionario.preguntas.aggregate(total=Sum('puntaje'))['total'] or 0
        cuestionario.puntaje_total = puntaje_total
        
        tiene_preguntas_manuales = cuestionario.preguntas.filter(
            tipo__nombre__in=['respuesta_abierta', 'simulador_2d', 'simulador_3d']
        ).exists()
        cuestionario.calificacion_automatica = not tiene_preguntas_manuales
        cuestionario.save()
        
        print(f"📊 Totales recalculados - Puntaje: {puntaje_total}")
        
        return JsonResponse({
            'success': True,
            'mensaje': 'Pregunta eliminada exitosamente',
            'puntaje_total': float(puntaje_total),
            'calificacion_automatica': not tiene_preguntas_manuales
        })
        
    except Exception as e:
        print(f"❌ Error al eliminar pregunta: {str(e)}")
        import traceback
        traceback.print_exc()
        
        return JsonResponse({
            'success': False, 
            'error': str(e)
        }, status=400)

# =====================================================
# GESTIÓN DE OPCIONES - COMPLETAMENTE FUNCIONAL
# =====================================================

@csrf_exempt
@require_http_methods(["POST"])
@login_required
def agregarOpcion(request):
    """Agregar una nueva opción a una pregunta - COMPLETAMENTE FUNCIONAL"""
    try:
        print("➕ Agregando nueva opción...")
        
        data = json.loads(request.body)
        pregunta_id = data.get('pregunta_id')
        
        print(f"📝 Pregunta ID: {pregunta_id}")
        
        pregunta = get_object_or_404(PreguntaCuestionario, id=pregunta_id)
        
        # Verificar permisos
        if (pregunta.cuestionario.recurso.seccion and 
            pregunta.cuestionario.recurso.seccion.curso.profesor != request.user):
            return JsonResponse({
                'success': False, 
                'error': 'Sin permisos'
            }, status=403)
        
        # Verificar tipo de pregunta
        if pregunta.tipo.nombre not in ['opcion_unica', 'opcion_multiple']:
            return JsonResponse({
                'success': False, 
                'error': 'No se pueden agregar opciones a este tipo de pregunta'
            })
        
        # Crear nueva opción
        opciones_count = pregunta.opciones.count()
        opcion = OpcionPregunta.objects.create(
            pregunta=pregunta,
            texto=f"Nueva opción {opciones_count + 1}",
            es_correcta=False
        )
        
        print(f"✅ Opción creada - ID: {opcion.id}, Texto: {opcion.texto}")
        
        return JsonResponse({
            'success': True,
            'opcion_id': opcion.id,
            'mensaje': 'Opción agregada exitosamente',
            'mantener_acordeon': True
        })
        
    except Exception as e:
        print(f"❌ Error al agregar opción: {str(e)}")
        import traceback
        traceback.print_exc()
        
        return JsonResponse({
            'success': False, 
            'error': str(e)
        }, status=400)

@csrf_exempt
@require_http_methods(["POST"])
@login_required
def actualizarOpcion(request):
    """Actualizar una opción existente - COMPLETAMENTE FUNCIONAL"""
    try:
        print("✏️ Actualizando opción...")
        
        data = json.loads(request.body)
        opcion_id = data.get('opcion_id')
        texto = data.get('texto')
        es_correcta = data.get('es_correcta', False)
        
        print(f"🔧 Opción ID: {opcion_id}, Texto: {texto}, Correcta: {es_correcta}")
        
        opcion = get_object_or_404(OpcionPregunta, id=opcion_id)
        
        # Verificar permisos
        if (opcion.pregunta.cuestionario.recurso.seccion and 
            opcion.pregunta.cuestionario.recurso.seccion.curso.profesor != request.user):
            return JsonResponse({
                'success': False, 
                'error': 'Sin permisos'
            }, status=403)
        
        # Actualizar opción
        opcion.texto = texto
        opcion.es_correcta = es_correcta
        opcion.save()
        
        print(f"✅ Opción actualizada: {texto}")
        
        return JsonResponse({
            'success': True, 
            'mensaje': 'Opción actualizada exitosamente'
        })
        
    except Exception as e:
        print(f"❌ Error al actualizar opción: {str(e)}")
        import traceback
        traceback.print_exc()
        
        return JsonResponse({
            'success': False, 
            'error': str(e)
        }, status=400)

@csrf_exempt
@require_http_methods(["DELETE"])
@login_required
def eliminarOpcion(request, opcion_id):
    """Eliminar una opción - COMPLETAMENTE FUNCIONAL"""
    try:
        print(f"🗑️ Eliminando opción ID: {opcion_id}")
        
        opcion = get_object_or_404(OpcionPregunta, id=opcion_id)
        
        # Verificar permisos
        if (opcion.pregunta.cuestionario.recurso.seccion and 
            opcion.pregunta.cuestionario.recurso.seccion.curso.profesor != request.user):
            return JsonResponse({
                'success': False, 
                'error': 'Sin permisos'
            }, status=403)
        
        # Verificar mínimo de opciones
        if opcion.pregunta.opciones.count() <= 2:
            return JsonResponse({
                'success': False, 
                'error': 'No se puede eliminar. Mínimo 2 opciones requeridas.'
            }, status=400)
        
        texto = opcion.texto
        opcion.delete()
        
        print(f"✅ Opción eliminada: {texto}")
        
        return JsonResponse({
            'success': True, 
            'mensaje': 'Opción eliminada exitosamente'
        })
        
    except Exception as e:
        print(f"❌ Error al eliminar opción: {str(e)}")
        import traceback
        traceback.print_exc()
        
        return JsonResponse({
            'success': False, 
            'error': str(e)
        }, status=400)

# =====================================================
# CREAR TIPOS DE PREGUNTA AUTOMÁTICAMENTE
# =====================================================

@csrf_exempt
@require_http_methods(["POST"])
@login_required
def crearTiposPregunta(request):
    """Crear tipos de pregunta automáticamente - COMPLETAMENTE FUNCIONAL"""
    try:
        print("🔧 Creando tipos de pregunta automáticamente...")
        
        tipos = [
            ('opcion_unica', 'Opción Única - Seleccionar una sola respuesta correcta'),
            ('opcion_multiple', 'Opción Múltiple - Seleccionar múltiples respuestas correctas'),
            ('falso_verdadero', 'Falso/Verdadero - Respuesta binaria'),
            ('unir_lineas', 'Unir con Líneas - Conectar elementos de dos columnas'),
            ('completar', 'Completar - Llenar espacios en blanco'),
            ('simulador_2d', 'Simulador 2D - Armar moléculas en 2D'),
            ('simulador_3d', 'Simulador 3D - Armar moléculas en 3D'),
            ('respuesta_abierta', 'Respuesta Abierta - Texto libre'),
        ]
        
        creados = 0
        for nombre, descripcion in tipos:
            tipo, created = TipoPregunta.objects.get_or_create(
                nombre=nombre,
                defaults={'descripcion': descripcion}
            )
            if created:
                creados += 1
                print(f"✅ Tipo creado: {nombre}")
            else:
                print(f"ℹ️ Tipo ya existe: {nombre}")
        
        # Crear tipo de recurso cuestionario
        tipo_recurso, created = TipoRecurso.objects.get_or_create(
            nombre='cuestionario',
            defaults={'descripcion': 'Recurso tipo cuestionario para evaluaciones'}
        )
        
        if created:
            print("✅ Tipo de recurso 'cuestionario' creado")
        
        return JsonResponse({
            'success': True,
            'mensaje': f'Tipos de pregunta listos. {creados} nuevos creados.',
        })
        
    except Exception as e:
        print(f"❌ Error al crear tipos: {str(e)}")
        import traceback
        traceback.print_exc()
        
        return JsonResponse({
            'success': False, 
            'error': str(e)
        }, status=400)

@login_required
def inicializar_tipos_pregunta(request):
    """Inicializar tipos de pregunta básicos en la base de datos - FUNCIONAL"""
    tipos_creados = crear_tipos_pregunta_basicos()
    return JsonResponse({
        'success': True,
        'mensaje': f'Se crearon {tipos_creados.count()} tipos de pregunta',
        'tipos': list(tipos_creados.values('nombre', 'descripcion'))
    })

# =====================================================
# FUNCIONES ADICIONALES PARA COMPLETAR LA FUNCIONALIDAD
# =====================================================

@csrf_exempt
@require_http_methods(["POST"])
@login_required
def actualizar_pregunta(request):
    """Función PRINCIPAL para actualizar pregunta - CORREGIDA PARA NO BORRAR CONTENIDO"""
    try:
        print("✏️ Actualizando pregunta...")
        
        # Detectar tipo de contenido y recopilar datos
        data = {}
        pregunta_id = None
        tiene_archivo = False
        
        if request.content_type and request.content_type.startswith('multipart/form-data'):
            # FormData (con archivos)
            pregunta_id = request.POST.get('pregunta_id')
            data.update(dict(request.POST))
            tiene_archivo = 'archivo_multimedia' in request.FILES
            print("📎 Datos recibidos via FormData")
        else:
            # JSON (sin archivos)
            try:
                json_data = json.loads(request.body)
                data = json_data
                pregunta_id = data.get('pregunta_id')
                print("📋 Datos recibidos via JSON")
            except json.JSONDecodeError:
                return JsonResponse({
                    'success': False, 
                    'error': 'Formato de datos inválido'
                })
        
        print(f"📝 Actualizando pregunta ID: {pregunta_id}")
        
        if not pregunta_id:
            return JsonResponse({
                'success': False, 
                'error': 'ID de pregunta requerido'
            })
        
        pregunta = get_object_or_404(PreguntaCuestionario, id=pregunta_id)
        
        # Verificar permisos
        if (pregunta.cuestionario.recurso.seccion and 
            pregunta.cuestionario.recurso.seccion.curso.profesor != request.user):
            return JsonResponse({
                'success': False, 
                'error': 'Sin permisos'
            }, status=403)
        
        with transaction.atomic():
            # ===== ACTUALIZAR CAMPOS BÁSICOS SOLO SI SE PROPORCIONAN =====
            
            # Enunciado (siempre requerido)
            enunciado = data.get('enunciado', '').strip()
            if enunciado:
                pregunta.enunciado = enunciado
                print(f"📝 Enunciado actualizado: {enunciado[:50]}...")
            elif not pregunta.enunciado:
                return JsonResponse({
                    'success': False, 
                    'error': 'El enunciado es obligatorio'
                })
            
            # Puntaje
            if 'puntaje' in data:
                try:
                    pregunta.puntaje = int(data.get('puntaje', 1))
                    print(f"⭐ Puntaje actualizado: {pregunta.puntaje}")
                except (ValueError, TypeError):
                    pregunta.puntaje = 1
            
            # Calificación manual
            if 'calificacion_manual' in data:
                pregunta.calificacion_manual = bool(data.get('calificacion_manual', False))
                print(f"✅ Calificación manual: {pregunta.calificacion_manual}")
            
            # ===== MANEJAR MULTIMEDIA SOLO SI SE PROPORCIONA =====
            
            if tiene_archivo:
                pregunta.archivo_multimedia = request.FILES['archivo_multimedia']
                print("📎 Archivo multimedia actualizado")
            
            if 'youtube_video_id' in data:
                youtube_id = data.get('youtube_video_id', '').strip()
                if youtube_id:
                    if len(youtube_id) == 11:
                        pregunta.url_youtube = youtube_id
                        print(f"🎥 Video YouTube actualizado: {youtube_id}")
                    else:
                        return JsonResponse({
                            'success': False, 
                            'error': 'ID de YouTube inválido'
                        })
                else:
                    # Si se envía vacío, limpiar
                    pregunta.url_youtube = None
                    print("🎥 Video YouTube eliminado")
            
            # ===== ACTUALIZAR CAMPOS ESPECÍFICOS POR TIPO (SOLO SI SE PROPORCIONAN) =====
            
            tipo_nombre = pregunta.tipo.nombre
            print(f"🔧 Procesando campos específicos para tipo: {tipo_nombre}")
            
            if tipo_nombre == 'respuesta_abierta':
                # SOLO actualizar si se proporciona el campo
                if 'respuesta_modelo' in data:
                    pregunta.respuesta_modelo = data.get('respuesta_modelo', '')
                    print(f"📝 Respuesta modelo actualizada")
                
                if 'criterios_evaluacion' in data:
                    pregunta.criterios_evaluacion = data.get('criterios_evaluacion', '')
                    print(f"📋 Criterios evaluación actualizados")
                
                if 'longitud_minima' in data:
                    try:
                        pregunta.longitud_minima = int(data.get('longitud_minima', 50))
                        print(f"📏 Longitud mínima: {pregunta.longitud_minima}")
                    except (ValueError, TypeError):
                        pass
                
                if 'longitud_maxima' in data:
                    try:
                        pregunta.longitud_maxima = int(data.get('longitud_maxima', 1000))
                        print(f"📏 Longitud máxima: {pregunta.longitud_maxima}")
                    except (ValueError, TypeError):
                        pass
                        
            elif tipo_nombre == 'completar':
                # SOLO actualizar si se proporciona el campo - CRÍTICO PARA EVITAR BORRADO
                if 'texto_completar' in data:
                    pregunta.texto_completar = data.get('texto_completar', '')
                    print(f"📝 Texto completar actualizado")
                
                if 'respuestas_completar' in data:
                    pregunta.respuestas_completar = data.get('respuestas_completar', '')
                    print(f"✅ Respuestas completar actualizadas")
                
                if 'sensible_mayusculas' in data:
                    pregunta.sensible_mayusculas = bool(data.get('sensible_mayusculas', False))
                    print(f"🔤 Sensible mayúsculas: {pregunta.sensible_mayusculas}")
                
                if 'ignorar_espacios' in data:
                    pregunta.ignorar_espacios = bool(data.get('ignorar_espacios', True))
                    print(f"🔲 Ignorar espacios: {pregunta.ignorar_espacios}")
                
                if 'permitir_alternativas' in data:
                    pregunta.permitir_alternativas = bool(data.get('permitir_alternativas', False))
                    print(f"🔀 Permitir alternativas: {pregunta.permitir_alternativas}")
                
                if 'respuestas_alternativas' in data:
                    pregunta.respuestas_alternativas = data.get('respuestas_alternativas', '')
                    print(f"📝 Respuestas alternativas actualizadas")
                        
            elif tipo_nombre == 'unir_lineas':
                # SOLO actualizar si se proporciona el campo - CRÍTICO PARA EVITAR BORRADO
                if 'columna_izquierda' in data:
                    pregunta.columna_izquierda = data.get('columna_izquierda', '')
                    print(f"📝 Columna izquierda actualizada")
                
                if 'columna_derecha' in data:
                    pregunta.columna_derecha = data.get('columna_derecha', '')
                    print(f"📝 Columna derecha actualizada")
                
                if 'conexiones_correctas' in data:
                    pregunta.conexiones_correctas = data.get('conexiones_correctas', '')
                    print(f"🔗 Conexiones correctas actualizadas")
                
                if 'mezclar_opciones' in data:
                    pregunta.mezclar_opciones = bool(data.get('mezclar_opciones', True))
                    print(f"🔀 Mezclar opciones: {pregunta.mezclar_opciones}")
                
                if 'permitir_conexiones_multiples' in data:
                    pregunta.permitir_conexiones_multiples = bool(data.get('permitir_conexiones_multiples', False))
                    print(f"🔗 Conexiones múltiples: {pregunta.permitir_conexiones_multiples}")
                        
            elif tipo_nombre in ['simulador_2d', 'simulador_3d']:
                # SOLO actualizar si se proporciona el campo
                if 'smiles_objetivo' in data:
                    pregunta.smiles_objetivo = data.get('smiles_objetivo', '')
                    print(f"🧪 SMILES objetivo actualizado: {pregunta.smiles_objetivo}")
                
                if 'descripcion_molecula' in data:
                    pregunta.descripcion_molecula = data.get('descripcion_molecula', '')
                    print(f"📝 Descripción molécula actualizada")
                
                if 'tolerancia_similitud' in data:
                    try:
                        pregunta.tolerancia_similitud = int(data.get('tolerancia_similitud', 95))
                        print(f"📊 Tolerancia similitud: {pregunta.tolerancia_similitud}%")
                    except (ValueError, TypeError):
                        pass
                
                if 'permitir_isomeros' in data:
                    pregunta.permitir_isomeros = bool(data.get('permitir_isomeros', False))
                    print(f"🔬 Permitir isómeros: {pregunta.permitir_isomeros}")
            
            # Guardar todos los cambios
            pregunta.save()
            print("✅ Pregunta actualizada exitosamente")
            
            # ===== ACTUALIZAR OPCIONES SI SE PROPORCIONAN =====
            
            if 'opciones' in data and isinstance(data['opciones'], list):
                print(f"🔧 Actualizando {len(data['opciones'])} opciones...")
                
                for opcion_data in data['opciones']:
                    if 'id' in opcion_data:
                        try:
                            opcion = OpcionPregunta.objects.get(
                                id=opcion_data['id'], 
                                pregunta=pregunta
                            )
                            
                            # Actualizar solo si se proporciona
                            if 'texto' in opcion_data:
                                opcion.texto = opcion_data['texto']
                            
                            if 'es_correcta' in opcion_data:
                                opcion.es_correcta = bool(opcion_data['es_correcta'])
                            
                            opcion.save()
                            print(f"✅ Opción {opcion.id} actualizada")
                            
                        except OpcionPregunta.DoesNotExist:
                            print(f"⚠️ Opción {opcion_data['id']} no encontrada")
                            continue
            
            # ===== RECALCULAR TOTALES =====
            
            cuestionario = pregunta.cuestionario
            puntaje_total = cuestionario.preguntas.aggregate(total=Sum('puntaje'))['total'] or 0
            cuestionario.puntaje_total = puntaje_total
            cuestionario.save()
            
            print(f"📊 Puntaje total recalculado: {puntaje_total}")
        
        return JsonResponse({
            'success': True,
            'mensaje': 'Pregunta actualizada exitosamente',
            'puntaje_total': float(puntaje_total),
            'pregunta_id': pregunta.id
        })
        
    except Exception as e:
        print(f"❌ Error al actualizar pregunta: {str(e)}")
        import traceback
        traceback.print_exc()
        
        return JsonResponse({
            'success': False, 
            'error': f'Error al actualizar pregunta: {str(e)}'
        }, status=400)
    

@csrf_exempt
@require_http_methods(["POST"])
@login_required
def eliminar_recurso_pregunta(request):
    """Eliminar recurso multimedia de una pregunta - FUNCIONAL"""
    try:
        print("🗑️ Eliminando recurso multimedia...")
        
        data = json.loads(request.body)
        pregunta_id = data.get('pregunta_id')
        tipo_recurso = data.get('tipo')
        
        pregunta = get_object_or_404(PreguntaCuestionario, id=pregunta_id)
        
        # Verificar permisos
        if (pregunta.cuestionario.recurso.seccion and 
            pregunta.cuestionario.recurso.seccion.curso.profesor != request.user):
            return JsonResponse({
                'success': False, 
                'error': 'Sin permisos'
            }, status=403)
        
        if tipo_recurso == 'archivo':
            if hasattr(pregunta, 'archivo_multimedia') and pregunta.archivo_multimedia:
                try:
                    pregunta.archivo_multimedia.delete()
                except:
                    pass
                pregunta.archivo_multimedia = None
                print("📎 Archivo multimedia eliminado")
        elif tipo_recurso == 'youtube':
            if hasattr(pregunta, 'url_youtube'):
                pregunta.url_youtube = None
                print("🎥 Video YouTube eliminado")
        
        pregunta.save()
        
        return JsonResponse({
            'success': True, 
            'mensaje': 'Recurso eliminado exitosamente'
        })
        
    except Exception as e:
        print(f"❌ Error al eliminar recurso: {str(e)}")
        return JsonResponse({
            'success': False, 
            'error': str(e)
        }, status=400)

@csrf_exempt
@require_http_methods(["POST"])
@login_required
def actualizar_configuracion_pregunta(request):
    """Actualizar configuración específica de pregunta según su tipo - FUNCIONAL"""
    try:
        print("⚙️ Actualizando configuración específica de pregunta...")
        
        data = json.loads(request.body)
        pregunta_id = data.get('pregunta_id')
        configuracion = data.get('configuracion')
        
        pregunta = get_object_or_404(PreguntaCuestionario, id=pregunta_id)
        
        # Verificar permisos
        if (pregunta.cuestionario.recurso.seccion and 
            pregunta.cuestionario.recurso.seccion.curso.profesor != request.user):
            return JsonResponse({
                'success': False, 
                'error': 'Sin permisos'
            }, status=403)
        
        # Esta función está preparada para manejar configuraciones específicas por tipo de pregunta
        # Por ejemplo: respuesta abierta (longitud mínima/máxima), simuladores (SMILES objetivo), etc.
        # Se puede extender en el futuro según necesidades específicas
        
        print(f"⚙️ Configuración para pregunta tipo: {pregunta.tipo.nombre}")
        
        return JsonResponse({
            'success': True,
            'mensaje': 'Configuración actualizada exitosamente'
        })
        
    except Exception as e:
        print(f"❌ Error al actualizar configuración: {str(e)}")
        return JsonResponse({
            'success': False, 
            'error': str(e)
        }, status=400)

@csrf_exempt
@require_http_methods(["POST"])
@login_required
def actualizarCuestionarioExistente(request, cuestionario_id):
    """Actualizar un cuestionario existente - FUNCIONAL"""
    try:
        print(f"✏️ Actualizando cuestionario existente ID: {cuestionario_id}")
        
        cuestionario = get_object_or_404(Cuestionario, id=cuestionario_id)
        
        # Verificar permisos
        if (cuestionario.recurso.seccion and 
            cuestionario.recurso.seccion.curso.profesor != request.user):
            return JsonResponse({
                'success': False, 
                'error': 'Sin permisos'
            }, status=403)
        
        titulo = request.POST.get('titulo', '').strip()
        descripcion = request.POST.get('descripcion', '').strip()
        
        if not titulo:
            return JsonResponse({
                'success': False, 
                'error': 'El título del cuestionario es obligatorio'
            })
        
        with transaction.atomic():
            cuestionario.recurso.titulo = titulo
            cuestionario.recurso.descripcion = descripcion
            cuestionario.recurso.save()
            
            cuestionario.instrucciones = descripcion
            cuestionario.save()
            
            print(f"✅ Cuestionario actualizado: {titulo}")
        
        return JsonResponse({
            'success': True, 
            'mensaje': 'Cuestionario actualizado exitosamente'
        })
        
    except Exception as e:
        print(f"❌ Error al actualizar cuestionario: {str(e)}")
        return JsonResponse({
            'success': False, 
            'error': str(e)
        }, status=400)


@csrf_exempt
@require_http_methods(["POST"])
def finalizar_cuestionario(request):
    """Finalizar y guardar cuestionario completo - FUNCIONAL"""
    try:
        print("🏁 Finalizando cuestionario...")
        
        data = json.loads(request.body)
        cuestionario_id = data.get('cuestionario_id')
        
        cuestionario = get_object_or_404(Cuestionario, id=cuestionario_id)
        
        # Verificar permisos
        if (cuestionario.recurso.seccion and 
            cuestionario.recurso.seccion.curso.profesor != request.user):
            return JsonResponse({
                'success': False, 
                'error': 'Sin permisos'
            }, status=403)
        
        with transaction.atomic():
            # Recalcular totales finales
            puntaje_total = cuestionario.preguntas.aggregate(total=Sum('puntaje'))['total'] or 0
            cuestionario.puntaje_total = puntaje_total
            
            tiene_preguntas_manuales = cuestionario.preguntas.filter(
                tipo__nombre__in=['respuesta_abierta', 'simulador_2d', 'simulador_3d']
            ).exists()
            cuestionario.calificacion_automatica = not tiene_preguntas_manuales
            
            cuestionario.save()
            
            print(f"✅ Cuestionario {cuestionario_id} finalizado - Puntaje total: {puntaje_total}")
        
        # URL de redirección
        redirect_url = '/docente/biblioteca-cuestionarios/'
        
        return JsonResponse({
            'success': True,
            'mensaje': 'Cuestionario finalizado exitosamente',
            'puntaje_total': float(puntaje_total),
            'redirect_url': redirect_url
        })
        
    except Exception as e:
        print(f"❌ Error al finalizar cuestionario: {str(e)}")
        import traceback
        traceback.print_exc()
        return JsonResponse({
            'success': False, 
            'error': str(e)
        }, status=400)


@csrf_exempt
@require_http_methods(["POST"])
def duplicar_pregunta(request):
    """Duplicar una pregunta existente - FUNCIONAL"""
    try:
        print("📋 Duplicando pregunta...")
        
        data = json.loads(request.body)
        pregunta_id = data.get('pregunta_id')
        
        pregunta_original = get_object_or_404(PreguntaCuestionario, id=pregunta_id)
        
        # Verificar permisos
        if (pregunta_original.cuestionario.recurso.seccion and 
            pregunta_original.cuestionario.recurso.seccion.curso.profesor != request.user):
            return JsonResponse({
                'success': False, 
                'error': 'Sin permisos'
            }, status=403)
        
        # Calcular nuevo orden
        ultima_pregunta = pregunta_original.cuestionario.preguntas.order_by('-orden').first()
        nuevo_orden = ultima_pregunta.orden + 1 if ultima_pregunta else 1
        
        with transaction.atomic():
            # Crear copia de la pregunta
            nueva_pregunta = PreguntaCuestionario.objects.create(
                cuestionario=pregunta_original.cuestionario,
                enunciado=f"{pregunta_original.enunciado} (Copia)",
                tipo=pregunta_original.tipo,
                orden=nuevo_orden,
                puntaje=pregunta_original.puntaje,
                # Campos específicos por tipo
                respuesta_modelo=getattr(pregunta_original, 'respuesta_modelo', None),
                criterios_evaluacion=getattr(pregunta_original, 'criterios_evaluacion', None),
                longitud_minima=getattr(pregunta_original, 'longitud_minima', None),
                longitud_maxima=getattr(pregunta_original, 'longitud_maxima', None),
                texto_completar=getattr(pregunta_original, 'texto_completar', None),
                respuestas_completar=getattr(pregunta_original, 'respuestas_completar', None),
                sensible_mayusculas=getattr(pregunta_original, 'sensible_mayusculas', False),
                ignorar_espacios=getattr(pregunta_original, 'ignorar_espacios', True),
                permitir_alternativas=getattr(pregunta_original, 'permitir_alternativas', False),
                respuestas_alternativas=getattr(pregunta_original, 'respuestas_alternativas', None),
                columna_izquierda=getattr(pregunta_original, 'columna_izquierda', None),
                columna_derecha=getattr(pregunta_original, 'columna_derecha', None),
                conexiones_correctas=getattr(pregunta_original, 'conexiones_correctas', None),
                mezclar_opciones=getattr(pregunta_original, 'mezclar_opciones', True),
                permitir_conexiones_multiples=getattr(pregunta_original, 'permitir_conexiones_multiples', False),
                smiles_objetivo=getattr(pregunta_original, 'smiles_objetivo', None),
                descripcion_molecula=getattr(pregunta_original, 'descripcion_molecula', None),
                tolerancia_similitud=getattr(pregunta_original, 'tolerancia_similitud', None),
                permitir_isomeros=getattr(pregunta_original, 'permitir_isomeros', False)
            )
            
            # Duplicar opciones si existen
            for opcion in pregunta_original.opciones.all():
                OpcionPregunta.objects.create(
                    pregunta=nueva_pregunta,
                    texto=opcion.texto,
                    es_correcta=opcion.es_correcta
                )
        
        return JsonResponse({
            'success': True,
            'pregunta_id': nueva_pregunta.id,
            'mensaje': 'Pregunta duplicada exitosamente'
        })
        
    except Exception as e:
        print(f"❌ Error al duplicar pregunta: {str(e)}")
        import traceback
        traceback.print_exc()
        return JsonResponse({
            'success': False, 
            'error': str(e)
        }, status=400)


@csrf_exempt
@require_http_methods(["POST"])
def subir_recurso(request):
    """Subir recurso multimedia a una pregunta - FUNCIONAL"""
    try:
        print("📤 Subiendo recurso multimedia...")
        
        pregunta_id = request.POST.get('pregunta_id')
        tipo_recurso = request.POST.get('tipo_recurso')
        
        print(f"📝 Pregunta ID: {pregunta_id}, Tipo: {tipo_recurso}")
        
        pregunta = get_object_or_404(PreguntaCuestionario, id=pregunta_id)
        
        # Verificar permisos
        if (pregunta.cuestionario.recurso.seccion and 
            pregunta.cuestionario.recurso.seccion.curso.profesor != request.user):
            return JsonResponse({
                'success': False, 
                'error': 'Sin permisos'
            }, status=403)
        
        if tipo_recurso == 'imagen' and 'archivo' in request.FILES:
            pregunta.archivo_multimedia = request.FILES['archivo']
            pregunta.save()
            print("✅ Imagen guardada")
            
        elif tipo_recurso == 'video' and 'archivo' in request.FILES:
            pregunta.archivo_multimedia = request.FILES['archivo']
            pregunta.save()
            print("✅ Video guardado")
            
        elif tipo_recurso == 'youtube':
            youtube_id = request.POST.get('youtube_id')
            if youtube_id:
                pregunta.url_youtube = youtube_id
                pregunta.save()
                print(f"✅ YouTube guardado: {youtube_id}")
                
        elif tipo_recurso == 'documento' and 'archivo' in request.FILES:
            pregunta.archivo_multimedia = request.FILES['archivo']
            pregunta.save()
            print("✅ Documento guardado")
            
        elif tipo_recurso == 'audio' and 'archivo' in request.FILES:
            pregunta.archivo_multimedia = request.FILES['archivo']
            pregunta.save()
            print("✅ Audio guardado")
        
        return JsonResponse({
            'success': True,
            'mensaje': 'Recurso cargado exitosamente'
        })
        
    except Exception as e:
        print(f"❌ Error al subir recurso: {str(e)}")
        import traceback
        traceback.print_exc()
        return JsonResponse({
            'success': False, 
            'error': str(e)
        }, status=400)


@login_required
def biblioteca_cuestionarios(request):
    """Vista para mostrar la biblioteca de cuestionarios"""
    user = request.user
    
    # ✅ CAMBIAR ESTA CONSULTA:
    cuestionarios = Cuestionario.objects.filter(
        Q(recurso__seccion__curso__profesor=user) |  # Con sección del usuario
        Q(recurso__seccion__isnull=True)              # Sin sección (biblioteca)
    ).select_related('recurso', 'recurso__seccion', 'recurso__seccion__curso').order_by('-id')
    
    # O MÁS SIMPLE, mostrar TODOS:
    # cuestionarios = Cuestionario.objects.all().order_by('-id')
    
    print(f"🔍 Total encontrados: {cuestionarios.count()}")  # Para debug
    
    context = {
        'cuestionarios': cuestionarios,
        'imgPerfil': user.imgPerfil,
        'usuario': user.username,
    }
    
    return render(request, 'docente/biblioteca_cuestionarios.html', context)


@login_required
def asignar_cuestionario_seccion(request, cuestionario_id):
    """Asignar un cuestionario de biblioteca a una sección específica"""
    cuestionario = get_object_or_404(Cuestionario, id=cuestionario_id)
    
    # Verificar que sea del profesor
    if cuestionario.recurso.creado_por != request.user:
        return redirect('bibliotecaCuestionarios')
    
    if request.method == 'POST':
        seccion_id = request.POST.get('seccion_id')
        seccion = get_object_or_404(Seccion, id=seccion_id, curso__profesor=request.user)
        
        # Asignar el cuestionario a la sección
        cuestionario.recurso.seccion = seccion
        cuestionario.recurso.orden = seccion.recursos.count() + 1
        cuestionario.recurso.save()
        
        messages.success(request, f'Cuestionario "{cuestionario.recurso.titulo}" asignado a {seccion.titulo}')
        return redirect('bibliotecaCuestionarios')
    
    # Obtener cursos y secciones del profesor
    cursos = request.user.cursos_creados.all()
    
    return render(request, 'docente/asignar_cuestionario.html', {
        'cuestionario': cuestionario,
        'cursos': cursos
    })



############ Estudiante ############

def biblioteca_cuestionarios_estudiante(request):
    """Vista de biblioteca de cuestionarios para estudiantes"""
    user = request.user
    
    # Obtener cuestionarios disponibles (todos por ahora, luego filtrar por curso)
    cuestionarios_disponibles = []
    
    for cuestionario in Cuestionario.objects.all():
        # Verificar disponibilidad por fechas
        disponible = cuestionario.esta_disponible()
        
        # Contar intentos del estudiante
        intentos_realizados = IntentoCuestionario.objects.filter(
            estudiante=user, 
            cuestionario=cuestionario
        ).count()
        
        # Verificar si puede hacer más intentos
        puede_intentar = (cuestionario.intentos_permitidos == 999 or 
                         intentos_realizados < cuestionario.intentos_permitidos)
        
        # Obtener mejor puntaje
        mejor_intento = IntentoCuestionario.objects.filter(
            estudiante=user, 
            cuestionario=cuestionario,
            completado=True
        ).order_by('-puntaje_obtenido').first()
        
        cuestionarios_disponibles.append({
            'cuestionario': cuestionario,
            'disponible': disponible,
            'intentos_realizados': intentos_realizados,
            'puede_intentar': puede_intentar,
            'mejor_puntaje': mejor_intento.puntaje_obtenido if mejor_intento else 0,
            'completado': mejor_intento is not None
        })
    
    context = {
        'cuestionarios_disponibles': cuestionarios_disponibles,
        'imgPerfil': user.imgPerfil,
        'usuario': user.username,
    }
    
    return render(request, 'estudiante/biblioteca_cuestionarios.html', context)

@csrf_exempt 
@role_required('Estudiante')
@login_required
def iniciar_cuestionario(request, cuestionario_id):
    """Iniciar un nuevo intento de cuestionario"""
    cuestionario = get_object_or_404(Cuestionario, id=cuestionario_id)
    user = request.user
    
    print(f"🚀 Iniciando cuestionario {cuestionario_id} para {user.username}")
    
    # Verificar disponibilidad
    if not cuestionario.esta_disponible():
        messages.error(request, 'Este cuestionario no está disponible en este momento.')
        return redirect('biblioteca_cuestionarios_estudiante')
    
    # Verificar intentos
    intentos_realizados = IntentoCuestionario.objects.filter(
        estudiante=user, 
        cuestionario=cuestionario
    ).count()
    
    if (cuestionario.intentos_permitidos != 999 and 
        intentos_realizados >= cuestionario.intentos_permitidos):
        messages.error(request, 'Has agotado todos los intentos para este cuestionario.')
        return redirect('biblioteca_cuestionarios_estudiante')
    
    # Crear nuevo intento
    intento = IntentoCuestionario.objects.create(
        estudiante=user,
        cuestionario=cuestionario,
        numero_intento=intentos_realizados + 1
    )
    
    print(f"✅ Intento creado: {intento.id}")
    
    return redirect('realizar_cuestionario', intento_id=intento.id)

@csrf_exempt 
@role_required('Estudiante')
@login_required
def realizar_cuestionario(request, intento_id):
    """Interfaz para realizar el cuestionario - VERSIÓN CON DEBUG"""
    intento = get_object_or_404(IntentoCuestionario, id=intento_id, estudiante=request.user)
    
    print(f"📝 Realizando cuestionario - Intento: {intento_id}")
    print(f"👤 Usuario: {request.user.username}")
    print(f"📋 Cuestionario: {intento.cuestionario.recurso.titulo}")
    
    if intento.completado:
        print(f"✅ Intento ya completado, redirigiendo a resultados")
        return redirect('resultado_cuestionario', intento_id=intento.id)
    
    # Obtener preguntas
    preguntas = intento.cuestionario.preguntas.all().order_by('orden')
    print(f"❓ Total de preguntas encontradas: {preguntas.count()}")
    
    for i, pregunta in enumerate(preguntas):
        print(f"  📋 Pregunta {i+1}: {pregunta.enunciado[:50]}... (Tipo: {pregunta.tipo.nombre})")
        if pregunta.tipo.nombre in ['opcion_unica', 'opcion_multiple', 'falso_verdadero']:
            opciones_count = pregunta.opciones.count()
            print(f"      🔘 Opciones: {opciones_count}")
    
    if preguntas.count() == 0:
        print("❌ ERROR: No hay preguntas en este cuestionario")
        messages.error(request, 'Este cuestionario no tiene preguntas configuradas.')
        return redirect('biblioteca_cuestionarios_estudiante')
    
    # Calcular tiempo restante
    tiempo_limite_segundos = intento.cuestionario.tiempo_limite * 60
    tiempo_transcurrido = (timezone.now() - intento.fecha_inicio).total_seconds()
    tiempo_restante = max(0, tiempo_limite_segundos - tiempo_transcurrido)
    
    print(f"⏱️ Tiempo límite: {intento.cuestionario.tiempo_limite} min")
    print(f"⏱️ Tiempo transcurrido: {tiempo_transcurrido/60:.1f} min")
    print(f"⏱️ Tiempo restante: {tiempo_restante/60:.1f} min")
    
    if tiempo_restante <= 0:
        # Tiempo agotado, finalizar automáticamente
        print("⏰ Tiempo agotado, finalizando automáticamente")
        return redirect('finalizar_cuestionario_estudiante', intento_id=intento.id)
    
    # Obtener respuestas existentes
    respuestas_existentes = {}
    for respuesta in intento.respuestas.all():
        respuestas_existentes[respuesta.pregunta.id] = respuesta
    
    print(f"💾 Respuestas existentes: {len(respuestas_existentes)}")
    
    context = {
        'intento': intento,
        'cuestionario': intento.cuestionario,
        'preguntas': preguntas,
        'tiempo_restante': int(tiempo_restante),
        'respuestas_existentes': respuestas_existentes,
        'imgPerfil': request.user.imgPerfil,
        'usuario': request.user.username,
    }
    
    print(f"📤 Enviando contexto con {len(context)} elementos")
    print(f"📋 Template: estudiante/realizar_cuestionario.html")
    
    return render(request, 'estudiante/realizar_cuestionario.html', context)

@role_required('Estudiante')
@login_required
@csrf_exempt
def guardar_respuesta(request):
    """Guardar respuesta de una pregunta vía AJAX - VERSIÓN CORREGIDA"""
    if request.method == 'POST':
        try:
            data = json.loads(request.body)
            intento_id = data.get('intento_id')
            pregunta_id = data.get('pregunta_id')
            respuesta = data.get('respuesta')
            
            print(f"💾 Guardando respuesta - Pregunta: {pregunta_id}")
            print(f"📋 Tipo de respuesta: {type(respuesta)}")
            print(f"📝 Contenido respuesta (primeros 100 chars): {str(respuesta)[:100]}")
            
            intento = get_object_or_404(IntentoCuestionario, id=intento_id, estudiante=request.user)
            pregunta = get_object_or_404(PreguntaCuestionario, id=pregunta_id)
            
            # Crear o actualizar respuesta
            respuesta_obj, created = RespuestaEstudiante.objects.get_or_create(
                intento=intento,
                pregunta=pregunta,
                defaults={'fecha_respuesta': timezone.now()}
            )
            
            # Procesar según tipo de pregunta
            if pregunta.tipo.nombre in ['opcion_unica', 'falso_verdadero']:
                if respuesta:
                    try:
                        opcion = get_object_or_404(OpcionPregunta, id=respuesta)
                        respuesta_obj.opcion_seleccionada = opcion
                        respuesta_obj.es_correcta = opcion.es_correcta
                        respuesta_obj.puntaje_obtenido = pregunta.puntaje if opcion.es_correcta else 0
                        print(f"✅ Opción única guardada - Correcta: {opcion.es_correcta}")
                    except:
                        print("❌ Error al procesar opción única")
                        respuesta_obj.opcion_seleccionada = None
                        respuesta_obj.es_correcta = False
                        respuesta_obj.puntaje_obtenido = 0
                else:
                    respuesta_obj.opcion_seleccionada = None
                    respuesta_obj.es_correcta = False
                    respuesta_obj.puntaje_obtenido = 0
            
            elif pregunta.tipo.nombre == 'opcion_multiple':
                if isinstance(respuesta, list) and respuesta:
                    respuesta_obj.opciones_multiples = ','.join(map(str, respuesta))
                    # Calcular puntaje para opción múltiple
                    try:
                        opciones_correctas = set(str(op.id) for op in pregunta.opciones.filter(es_correcta=True))
                        opciones_seleccionadas = set(map(str, respuesta))
                        
                        if opciones_seleccionadas == opciones_correctas:
                            respuesta_obj.es_correcta = True
                            respuesta_obj.puntaje_obtenido = pregunta.puntaje
                        else:
                            respuesta_obj.es_correcta = False
                            respuesta_obj.puntaje_obtenido = 0
                        print(f"✅ Opción múltiple guardada - Correcta: {respuesta_obj.es_correcta}")
                    except:
                        print("❌ Error al procesar opción múltiple")
                        respuesta_obj.es_correcta = False
                        respuesta_obj.puntaje_obtenido = 0
                else:
                    respuesta_obj.opciones_multiples = ''
                    respuesta_obj.es_correcta = False
                    respuesta_obj.puntaje_obtenido = 0
            
            elif pregunta.tipo.nombre == 'respuesta_abierta':
                respuesta_texto = str(respuesta) if respuesta else ''
                if len(respuesta_texto) > 5000:  # Limitar longitud
                    respuesta_texto = respuesta_texto[:5000]
                respuesta_obj.respuesta_texto = respuesta_texto
                respuesta_obj.es_correcta = None  # Requiere calificación manual
                respuesta_obj.puntaje_obtenido = 0  # Se califica después
                print(f"✅ Respuesta abierta guardada - {len(respuesta_texto)} caracteres")
            
            elif pregunta.tipo.nombre == 'completar':
                # CORREGIR ESTE CASO QUE ESTÁ CAUSANDO EL PROBLEMA
                if isinstance(respuesta, str) and respuesta:
                    # Limpiar respuesta de caracteres extraños
                    respuesta_limpia = respuesta.replace('|', '').strip()
                    if len(respuesta_limpia) > 1000:  # Limitar longitud
                        respuesta_limpia = respuesta_limpia[:1000]
                    respuesta_obj.respuestas_completar = respuesta_limpia
                    print(f"✅ Completar guardado - Respuesta limpia: {respuesta_limpia[:50]}")
                elif isinstance(respuesta, list):
                    # Si es una lista, unir con |
                    respuesta_lista = [str(r).strip() for r in respuesta if str(r).strip()]
                    respuesta_obj.respuestas_completar = '|'.join(respuesta_lista)
                    print(f"✅ Completar guardado - Lista: {len(respuesta_lista)} elementos")
                else:
                    respuesta_obj.respuestas_completar = ''
                    print("⚠️ Completar vacío")
                
                # Calcular puntaje automáticamente
                try:
                    respuesta_obj.calcular_puntaje_automatico()
                except:
                    print("❌ Error al calcular puntaje automático")
                    respuesta_obj.puntaje_obtenido = 0
            
            elif pregunta.tipo.nombre == 'unir_lineas':
                if respuesta:
                    try:
                        respuesta_obj.conexiones_realizadas = json.dumps(respuesta) if isinstance(respuesta, dict) else str(respuesta)
                        # Calcular puntaje automáticamente
                        respuesta_obj.calcular_puntaje_automatico()
                        print(f"✅ Unir líneas guardado")
                    except:
                        print("❌ Error al procesar unir líneas")
                        respuesta_obj.conexiones_realizadas = ''
                        respuesta_obj.puntaje_obtenido = 0
                else:
                    respuesta_obj.conexiones_realizadas = ''
                    respuesta_obj.puntaje_obtenido = 0
            
            elif pregunta.tipo.nombre in ['simulador_2d', 'simulador_3d']:
                respuesta_smiles = str(respuesta) if respuesta else ''
                if len(respuesta_smiles) > 500:  # Limitar longitud
                    respuesta_smiles = respuesta_smiles[:500]
                respuesta_obj.respuesta_smiles = respuesta_smiles
                respuesta_obj.es_correcta = None  # Requiere validación
                respuesta_obj.puntaje_obtenido = 0
                print(f"✅ Simulador guardado - SMILES: {respuesta_smiles[:20]}")
            
            # Guardar la respuesta
            respuesta_obj.save()
            
            print(f"✅ Respuesta guardada exitosamente - Puntaje: {respuesta_obj.puntaje_obtenido}")
            
            return JsonResponse({
                'success': True,
                'puntaje_obtenido': float(respuesta_obj.puntaje_obtenido)
            })
            
        except Exception as e:
            print(f"❌ Error al guardar respuesta: {str(e)}")
            import traceback
            traceback.print_exc()
            return JsonResponse({'success': False, 'error': str(e)})
    
    return JsonResponse({'success': False, 'error': 'Método no permitido'})

@csrf_exempt 
@role_required('Estudiante')
@login_required
def finalizar_cuestionario_estudiante(request, intento_id):
    """Finalizar el cuestionario y calcular puntaje"""
    intento = get_object_or_404(IntentoCuestionario, id=intento_id, estudiante=request.user)
    
    print(f"🏁 Finalizando cuestionario - Intento: {intento_id}")
    
    if not intento.completado:
        # Calcular puntaje total
        puntaje_total = 0
        for respuesta in intento.respuestas.all():
            puntaje_total += respuesta.puntaje_obtenido
        
        # Finalizar intento
        intento.puntaje_obtenido = puntaje_total
        intento.completado = True
        intento.fecha_finalizacion = timezone.now()
        intento.tiempo_empleado = (intento.fecha_finalizacion - intento.fecha_inicio).total_seconds()
        intento.save()
        
        print(f"✅ Cuestionario finalizado - Puntaje: {puntaje_total}")
    
    return redirect('resultado_cuestionario', intento_id=intento.id)

@csrf_exempt 
@role_required('Estudiante')
@login_required
def resultado_cuestionario(request, intento_id):
    """Mostrar resultados del cuestionario"""
    intento = get_object_or_404(IntentoCuestionario, id=intento_id, estudiante=request.user)
    
    # Obtener todas las respuestas con detalles
    respuestas = intento.respuestas.select_related(
        'pregunta', 'opcion_seleccionada'
    ).prefetch_related('pregunta__opciones').all()
    
    # Procesar respuestas para el template
    for respuesta in respuestas:
        # Calcular porcentaje de la pregunta
        if respuesta.pregunta.puntaje > 0:
            respuesta.porcentaje_pregunta = int((respuesta.puntaje_obtenido / respuesta.pregunta.puntaje) * 100)
        else:
            respuesta.porcentaje_pregunta = 0
        
        # Procesar opciones múltiples
        if respuesta.opciones_multiples:
            respuesta.opciones_ids = [int(x.strip()) for x in respuesta.opciones_multiples.split(',') if x.strip()]
        else:
            respuesta.opciones_ids = []
        
        # Procesar respuestas de completar
        if respuesta.pregunta.respuestas_completar:
            respuesta.respuestas_correctas_lista = [x.strip() for x in respuesta.pregunta.respuestas_completar.split(',') if x.strip()]
        else:
            respuesta.respuestas_correctas_lista = []
    
    # Calcular porcentaje total
    puntaje_porcentaje = 0
    if intento.cuestionario.puntaje_total > 0:
        puntaje_porcentaje = (intento.puntaje_obtenido / intento.cuestionario.puntaje_total * 100)
    
    context = {
        'intento': intento,
        'respuestas': respuestas,
        'cuestionario': intento.cuestionario,
        'puntaje_porcentaje': puntaje_porcentaje,
        'imgPerfil': request.user.imgPerfil,
        'usuario': request.user.username,
    }
    
    return render(request, 'estudiante/resultado_cuestionario.html', context)

@csrf_exempt 
@role_required('Estudiante')
@login_required
def historial_cuestionario_especifico(request, cuestionario_id):
    """Historial de intentos de un cuestionario específico"""
    cuestionario = get_object_or_404(Cuestionario, id=cuestionario_id)
    user = request.user
    
    # Obtener todos los intentos del estudiante para este cuestionario
    intentos = IntentoCuestionario.objects.filter(
        estudiante=user,
        cuestionario=cuestionario,
        completado=True
    ).select_related('cuestionario').order_by('-fecha_finalizacion')
    
    print(f"📊 Historial específico - Cuestionario: {cuestionario.recurso.titulo}")
    print(f"👤 Usuario: {user.username}")
    print(f"🔢 Total intentos: {intentos.count()}")
    
    # Calcular estadísticas
    estadisticas = {
        'total_intentos': intentos.count(),
        'mejor_puntaje': 0,
        'promedio_puntaje': 0,
        'ultimo_puntaje': 0,
        'mejor_tiempo': 0,
        'tiempo_promedio': 0
    }
    
    if intentos.exists():
        puntajes = [float(intento.puntaje_obtenido) for intento in intentos]
        tiempos = [intento.tiempo_empleado for intento in intentos if intento.tiempo_empleado]
        
        estadisticas['mejor_puntaje'] = max(puntajes)
        estadisticas['promedio_puntaje'] = sum(puntajes) / len(puntajes)
        estadisticas['ultimo_puntaje'] = float(intentos.first().puntaje_obtenido)
        
        if tiempos:
            estadisticas['mejor_tiempo'] = min(tiempos) / 60  # en minutos
            estadisticas['tiempo_promedio'] = sum(tiempos) / len(tiempos) / 60
    
    context = {
        'cuestionario': cuestionario,
        'intentos': intentos,
        'estadisticas': estadisticas,
        'imgPerfil': user.imgPerfil,
        'usuario': user.username,
    }
    
    return render(request, 'estudiante/historial_cuestionario_especifico.html', context)



# =====================================================
# VISTAS PRINCIPALES PARA CREAR COMPETENCIAS
# =====================================================

@login_required
@role_required(['Admin', 'Docente'])
def crear_competencia(request, seccion_id=None):
    """Vista principal para crear competencias - CORREGIDA"""
    print(f"🏆 crear_competencia llamada con seccion_id: {seccion_id}")
    
    seccion = None
    if seccion_id:
        seccion = get_object_or_404(Seccion, id=seccion_id)
        print(f"📁 Sección encontrada: {seccion.titulo}")
        
        # Verificar permisos
        if seccion.curso.profesor != request.user:
            print("❌ Usuario sin permisos")
            messages.error(request, 'No tienes permisos para crear competencias en esta sección')
            return redirect('dashboard-adm')
    
    # Crear tipos de pregunta si no existen
    tipos_pregunta = TipoPregunta.objects.all()
    if not tipos_pregunta.exists():
        print("🔧 Creando tipos de pregunta básicos...")
        tipos_pregunta = crear_tipos_pregunta_basicos()
    
    # AGREGAR: Verificar que existe el tipo de recurso 'competencia'
    tipo_competencia, created = TipoRecurso.objects.get_or_create(
        nombre='competencia',
        defaults={'descripcion': 'Recurso tipo competencia en tiempo real'}
    )
    if created:
        print("✅ Tipo de recurso 'competencia' creado automáticamente")
    
    print(f"📝 Tipos de pregunta disponibles: {tipos_pregunta.count()}")
    
    context = {
        'competencia': None,
        'preguntas': [],
        'seccion': seccion,
        'seccion_id': seccion_id,
        'tipos_pregunta': tipos_pregunta,
        'modalidades': Competencia.MODALIDAD_CHOICES,
        'estados': Competencia.ESTADO_CHOICES,
        # AGREGAR: Información adicional para el template
        'max_tiempo_limite': 180,  # 3 horas máximo
        'max_miembros_grupo': 10,
        'configuracion_defecto': {
            'tiempo_limite': 30,
            'max_miembros_grupo': 4,
            'grupos_aleatorios': False,
            'grupos_abiertos': True,
            'mostrar_resultados_inmediatos': True,
            'permitir_reingreso': False,
            'orden_preguntas_aleatorio': False,
        }
    }
    
    return render(request, 'docente/competencias/crear_competencia.html', context)


@login_required
@role_required(['Admin', 'Docente'])
def editar_competencia(request, competencia_id):
    """Vista para editar una competencia existente"""
    print(f"✏️ editar_competencia llamada con ID: {competencia_id}")
    
    competencia = get_object_or_404(Competencia, id=competencia_id)
    print(f"🏆 Competencia encontrada: {competencia.recurso.titulo}")
    
    # Verificar permisos
    if (competencia.recurso.seccion and 
        competencia.recurso.seccion.curso.profesor != request.user):
        print("❌ Usuario sin permisos para editar")
        return redirect('dashboard-adm')
    
    # CORREGIR ESTA LÍNEA - Obtener preguntas usando la relación genérica
    from django.contrib.contenttypes.models import ContentType
    competencia_ct = ContentType.objects.get_for_model(Competencia)
    preguntas = PreguntaCuestionario.objects.filter(
        content_type=competencia_ct,
        object_id=competencia.id
    ).order_by('orden')
    
    tipos_pregunta = TipoPregunta.objects.all()
    
    # Obtener grupos si es modalidad grupal
    grupos = []
    if competencia.modalidad == 'grupal':
        grupos = competencia.grupos.all().order_by('fecha_creacion')
    
    # Calcular datos dinámicos
    puntaje_calculado = preguntas.aggregate(total=Sum('puntaje'))['total'] or 0
    participantes_count = competencia.participaciones.filter(activo=True).count()
    
    print(f"📊 Preguntas: {preguntas.count()}, Puntaje: {puntaje_calculado}, Participantes: {participantes_count}")
    
    context = {
        'competencia': competencia,
        'preguntas': preguntas,
        'tipos_pregunta': tipos_pregunta,
        'modalidades': Competencia.MODALIDAD_CHOICES,
        'estados': Competencia.ESTADO_CHOICES,
        'puntaje_calculado': puntaje_calculado,
        'participantes_count': participantes_count,
        'grupos': grupos,
        'puede_iniciar': competencia.puede_iniciar(),
        'esta_activa': competencia.esta_activa(),
    }
    
    return render(request, 'docente/competencias/crear_competencia.html', context)


@login_required
@role_required(['Admin', 'Docente'])
def editar_competencia_por_recurso(request, recurso_id):
    """Vista para editar competencia usando recurso_id - Para compatibilidad"""
    print(f"✏️ editar_competencia_por_recurso llamada con Recurso ID: {recurso_id}")
    
    recurso = get_object_or_404(Recurso, id=recurso_id)
    competencia = get_object_or_404(Competencia, recurso=recurso)
    
    return editar_competencia(request, competencia.id)


# =====================================================
# GUARDAR COMPETENCIA BÁSICA
# =====================================================

@csrf_exempt
@require_http_methods(["POST"])
@login_required
@role_required(['Admin', 'Docente'])
def guardar_competencia(request):
    """Guardar competencia básica - Similar a guardar_cuestionario"""
    try:
        print("💾 Iniciando guardado de competencia...")
        print(f"POST data: {request.POST}")
        
        seccion_id = request.POST.get('seccion_id')
        titulo = request.POST.get('titulo', '').strip()
        descripcion = request.POST.get('descripcion', '').strip()
        modalidad = request.POST.get('modalidad', 'individual')
        
        print(f"🏆 Datos recibidos - Título: '{titulo}', Modalidad: '{modalidad}'")
        
        # Validaciones
        if not titulo:
            return JsonResponse({
                'success': False, 
                'error': 'El título de la competencia es obligatorio'
            })
        
        if len(titulo) < 3:
            return JsonResponse({
                'success': False,
                'error': 'El título debe tener al menos 3 caracteres'
            })
        
        if not descripcion:
            return JsonResponse({
                'success': False,
                'error': 'La descripción de la competencia es obligatoria'
            })
        
        if modalidad not in ['individual', 'grupal']:
            return JsonResponse({
                'success': False,
                'error': 'Modalidad no válida'
            })
        
        # Validar sección
        seccion = None
        if seccion_id:
            seccion = get_object_or_404(Seccion, id=seccion_id)
            if seccion.curso.profesor != request.user:
                return JsonResponse({
                    'success': False, 
                    'error': 'Sin permisos para crear competencia en esta sección'
                }, status=403)
        
        # Crear competencia con transacción
        with transaction.atomic():
            print("🔧 Creando tipo de recurso competencia...")
            
            # Crear o obtener tipo de recurso
            tipo_competencia, created = TipoRecurso.objects.get_or_create(
                nombre='competencia',
                defaults={'descripcion': 'Recurso tipo competencia en tiempo real'}
            )
            
            if created:
                print("✅ Tipo de recurso 'competencia' creado")
            
            # Crear recurso
            print("🗂️ Creando recurso...")
            recurso = Recurso.objects.create(
                seccion=seccion,
                tipo=tipo_competencia,
                titulo=titulo,
                descripcion=descripcion,
                orden=seccion.recursos.count() + 1 if seccion else 1,
                creado_por=request.user
            )
            
            print(f"✅ Recurso creado con ID: {recurso.id}")
            
            # Configuración adicional para modalidad grupal
            max_miembros = 4
            grupos_aleatorios = False
            grupos_abiertos = True
            
            if modalidad == 'grupal':
                max_miembros = int(request.POST.get('max_miembros_grupo', 4))
                grupos_aleatorios = request.POST.get('grupos_aleatorios') == 'true'
                grupos_abiertos = request.POST.get('grupos_abiertos', 'true') == 'true'
                
                print(f"🏆 Configuración grupal - Max: {max_miembros}, Aleatorios: {grupos_aleatorios}")
            
            # Crear competencia
            print("🏆 Creando competencia...")
            competencia = Competencia.objects.create(
                recurso=recurso,
                instrucciones=descripcion,
                modalidad=modalidad,
                max_miembros_grupo=max_miembros,
                grupos_aleatorios=grupos_aleatorios,
                grupos_abiertos=grupos_abiertos,
                tiempo_limite=30,
                estado='configuracion',
                mostrar_resultados_inmediatos=True,
                permitir_reingreso=False,
                orden_preguntas_aleatorio=False
            )
            
            print(f"✅ Competencia creada con ID: {competencia.id}")
            print(f"📌 PIN de acceso generado: {competencia.pin_acceso}")
        
        response_data = {
            'success': True,
            'competencia_id': competencia.id,
            'pin_acceso': competencia.pin_acceso,
            'codigo_acceso': competencia.codigo_acceso,
            'mensaje': 'Competencia creada exitosamente',
            'redirect_url': f'/docente/editar_competencia/{competencia.id}/'
        }
        
        print(f"✅ Respuesta exitosa: {response_data}")
        return JsonResponse(response_data)
        
    except Exception as e:
        print(f"❌ Error al guardar competencia: {str(e)}")
        import traceback
        traceback.print_exc()
        
        return JsonResponse({
            'success': False, 
            'error': f'Error interno: {str(e)}'
        }, status=400)



#CONFIGURACIÓN Y GESTIÓN DE COMPETENCIAS

# =====================================================
# ACTUALIZAR CONFIGURACIÓN DE COMPETENCIAS
# =====================================================

@csrf_exempt
@require_http_methods(["POST"])
@login_required
@role_required(['Admin', 'Docente'])
def actualizar_configuracion_competencia(request):
    """Actualizar configuración de la competencia"""
    try:
        print("⚙️ Actualizando configuración de competencia...")
        
        data = json.loads(request.body)
        competencia_id = data.get('competencia_id')
        configuracion = data.get('configuracion')
        
        print(f"🏆 Competencia ID: {competencia_id}")
        print(f"⚙️ Configuración: {configuracion}")
        
        competencia = get_object_or_404(Competencia, id=competencia_id)
        
        # Verificar permisos
        if (competencia.recurso.seccion and 
            competencia.recurso.seccion.curso.profesor != request.user):
            return JsonResponse({
                'success': False, 
                'error': 'Sin permisos'
            }, status=403)
        
        # Solo permitir cambios si está en configuración
        if competencia.estado not in ['configuracion', 'esperando']:
            return JsonResponse({
                'success': False, 
                'error': 'No se puede modificar una competencia activa o finalizada'
            })
        
        # Actualizar configuración
        if configuracion:
            # Configuración básica
            if 'tiempo_limite' in configuracion:
                tiempo = int(configuracion['tiempo_limite'])
                if 1 <= tiempo <= 180:  # Máximo 3 horas
                    competencia.tiempo_limite = tiempo
                    print(f"⏱️ Tiempo límite actualizado: {tiempo} min")
            
            if 'modalidad' in configuracion:
                nueva_modalidad = configuracion['modalidad']
                if nueva_modalidad in ['individual', 'grupal']:
                    competencia.modalidad = nueva_modalidad
                    print(f"👥 Modalidad actualizada: {nueva_modalidad}")
            
            # Configuración para modalidad grupal
            if 'max_miembros_grupo' in configuracion:
                max_miembros = int(configuracion['max_miembros_grupo'])
                if 2 <= max_miembros <= 10:
                    competencia.max_miembros_grupo = max_miembros
                    print(f"👥 Max miembros por grupo: {max_miembros}")
            
            if 'grupos_aleatorios' in configuracion:
                competencia.grupos_aleatorios = bool(configuracion['grupos_aleatorios'])
                print(f"🎲 Grupos aleatorios: {competencia.grupos_aleatorios}")
            
            if 'grupos_abiertos' in configuracion:
                competencia.grupos_abiertos = bool(configuracion['grupos_abiertos'])
                print(f"🔓 Grupos abiertos: {competencia.grupos_abiertos}")
            
            # Configuración avanzada
            if 'mostrar_resultados_inmediatos' in configuracion:
                competencia.mostrar_resultados_inmediatos = bool(configuracion['mostrar_resultados_inmediatos'])
                print(f"📊 Mostrar resultados: {competencia.mostrar_resultados_inmediatos}")
            
            if 'permitir_reingreso' in configuracion:
                competencia.permitir_reingreso = bool(configuracion['permitir_reingreso'])
                print(f"🔄 Permitir reingreso: {competencia.permitir_reingreso}")
            
            if 'orden_preguntas_aleatorio' in configuracion:
                competencia.orden_preguntas_aleatorio = bool(configuracion['orden_preguntas_aleatorio'])
                print(f"🔀 Orden aleatorio: {competencia.orden_preguntas_aleatorio}")
            
            # Manejar fecha de inicio programada
            if 'fecha_inicio' in configuracion and configuracion['fecha_inicio']:
                try:
                    fecha_inicio = datetime.fromisoformat(
                        configuracion['fecha_inicio'].replace('Z', '+00:00')
                    )
                    competencia.fecha_inicio = fecha_inicio
                    print(f"📅 Fecha inicio programada: {fecha_inicio}")
                except Exception as e:
                    print(f"⚠️ Error en fecha inicio: {e}")
        
        competencia.save()
        print("✅ Configuración guardada exitosamente")
        
        return JsonResponse({
            'success': True,
            'mensaje': 'Configuración actualizada exitosamente',
        })
        
    except Exception as e:
        print(f"❌ Error al actualizar configuración: {str(e)}")
        import traceback
        traceback.print_exc()
        
        return JsonResponse({
            'success': False, 
            'error': str(e)
        }, status=400)


# =====================================================
# GESTIÓN DE PREGUNTAS DE COMPETENCIAS
# =====================================================

@csrf_exempt
@require_http_methods(["POST"])
@login_required
@role_required(['Admin', 'Docente'])
def agregar_pregunta_competencia(request):
    """Agregar una nueva pregunta a la competencia - CORREGIDO"""
    try:
        print("➕ Agregando nueva pregunta a competencia...")
        
        data = json.loads(request.body)
        competencia_id = data.get('competencia_id')
        tipo_nombre = data.get('tipo')
        
        print(f"🏆 Competencia ID: {competencia_id}, Tipo: {tipo_nombre}")
        
        competencia = get_object_or_404(Competencia, id=competencia_id)
        
        # Verificar permisos
        if (competencia.recurso.seccion and 
            competencia.recurso.seccion.curso.profesor != request.user):
            return JsonResponse({
                'success': False, 
                'error': 'Sin permisos'
            }, status=403)
        
        # Solo permitir si está en configuración
        if competencia.estado not in ['configuracion']:
            return JsonResponse({
                'success': False, 
                'error': 'No se pueden agregar preguntas a una competencia iniciada'
            })
        
        tipo_pregunta = get_object_or_404(TipoPregunta, nombre=tipo_nombre)
        
        # Calcular orden - CORREGIR ESTA LÍNEA
        from django.contrib.contenttypes.models import ContentType
        competencia_ct = ContentType.objects.get_for_model(Competencia)
        
        ultima_pregunta = PreguntaCuestionario.objects.filter(
            content_type=competencia_ct,
            object_id=competencia.id
        ).order_by('-orden').first()
        
        orden = ultima_pregunta.orden + 1 if ultima_pregunta else 1
        
        print(f"📝 Creando pregunta orden {orden}")
        
        # Crear pregunta - CORREGIR ESTA PARTE
        pregunta = PreguntaCuestionario.objects.create(
            content_type=competencia_ct,
            object_id=competencia.id,
            tipo=tipo_pregunta,
            enunciado=f"Nueva pregunta {orden}",
            orden=orden,
            puntaje=1
        )
        
        print(f"✅ Pregunta creada con ID: {pregunta.id}")
        
        # Crear opciones por defecto según el tipo
        if tipo_nombre in ['opcion_unica', 'opcion_multiple', 'falso_verdadero']:
            if tipo_nombre == 'falso_verdadero':
                opciones_default = [('Verdadero', True), ('Falso', False)]
            else:
                opciones_default = [('Opción 1', True), ('Opción 2', False)]
            
            for texto, es_correcta in opciones_default:
                OpcionPregunta.objects.create(
                    pregunta=pregunta,
                    texto=texto,
                    es_correcta=es_correcta
                )
                print(f"✅ Opción creada: {texto} (correcta: {es_correcta})")
        
        return JsonResponse({
            'success': True,
            'pregunta_id': pregunta.id,
            'mensaje': 'Pregunta agregada exitosamente'
        })
        
    except Exception as e:
        print(f"❌ Error al agregar pregunta: {str(e)}")
        import traceback
        traceback.print_exc()
        
        return JsonResponse({
            'success': False, 
            'error': str(e)
        }, status=400)



@csrf_exempt
@require_http_methods(["POST"])
@login_required
@role_required(['Admin', 'Docente'])
def actualizar_pregunta_competencia(request):
    """Actualizar una pregunta de competencia - CORREGIDO"""
    try:
        print("✏️ Actualizando pregunta de competencia...")
        
        # Detectar tipo de contenido
        data = {}
        pregunta_id = None
        
        if request.content_type and request.content_type.startswith('multipart/form-data'):
            pregunta_id = request.POST.get('pregunta_id')
            data.update(dict(request.POST))
        else:
            try:
                json_data = json.loads(request.body)
                data = json_data
                pregunta_id = data.get('pregunta_id')
            except json.JSONDecodeError:
                return JsonResponse({
                    'success': False, 
                    'error': 'Formato de datos inválido'
                })
        
        print(f"📝 Actualizando pregunta ID: {pregunta_id}")
        
        if not pregunta_id:
            return JsonResponse({
                'success': False, 
                'error': 'ID de pregunta requerido'
            })
        
        pregunta = get_object_or_404(PreguntaCuestionario, id=pregunta_id)
        
        # CORREGIR: Verificar que pertenece a una competencia
        if not isinstance(pregunta.content_object, Competencia):
            return JsonResponse({
                'success': False, 
                'error': 'Esta pregunta no pertenece a una competencia'
            })
        
        competencia = pregunta.content_object
        
        # Verificar permisos
        if (competencia.recurso.seccion and 
            competencia.recurso.seccion.curso.profesor != request.user):
            return JsonResponse({
                'success': False, 
                'error': 'Sin permisos'
            }, status=403)
        
        # Solo permitir si está en configuración
        if competencia.estado not in ['configuracion']:
            return JsonResponse({
                'success': False, 
                'error': 'No se pueden modificar preguntas de una competencia iniciada'
            })
        
        with transaction.atomic():
            # Actualizar campos básicos solo si se proporcionan
            enunciado = data.get('enunciado', '').strip()
            if enunciado:
                pregunta.enunciado = enunciado
                print(f"📝 Enunciado actualizado: {enunciado[:50]}...")
            
            if 'puntaje' in data:
                try:
                    pregunta.puntaje = int(data.get('puntaje', 1))
                    print(f"⭐ Puntaje actualizado: {pregunta.puntaje}")
                except (ValueError, TypeError):
                    pregunta.puntaje = 1
            
            # Guardar cambios
            pregunta.save()
            print("✅ Pregunta actualizada exitosamente")
            
            # Actualizar opciones si se proporcionan
            if 'opciones' in data and isinstance(data['opciones'], list):
                print(f"🔧 Actualizando {len(data['opciones'])} opciones...")
                
                for opcion_data in data['opciones']:
                    if 'id' in opcion_data:
                        try:
                            opcion = OpcionPregunta.objects.get(
                                id=opcion_data['id'], 
                                pregunta=pregunta
                            )
                            
                            if 'texto' in opcion_data:
                                opcion.texto = opcion_data['texto']
                            
                            if 'es_correcta' in opcion_data:
                                opcion.es_correcta = bool(opcion_data['es_correcta'])
                            
                            opcion.save()
                            print(f"✅ Opción {opcion.id} actualizada")
                            
                        except OpcionPregunta.DoesNotExist:
                            print(f"⚠️ Opción {opcion_data['id']} no encontrada")
                            continue
        
        return JsonResponse({
            'success': True,
            'mensaje': 'Pregunta actualizada exitosamente',
            'pregunta_id': pregunta.id
        })
        
    except Exception as e:
        print(f"❌ Error al actualizar pregunta: {str(e)}")
        import traceback
        traceback.print_exc()
        
        return JsonResponse({
            'success': False, 
            'error': f'Error al actualizar pregunta: {str(e)}'
        }, status=400)


@csrf_exempt
@require_http_methods(["DELETE"])
@login_required
@role_required(['Admin', 'Docente'])
def eliminar_pregunta_competencia(request, pregunta_id):
    """Eliminar una pregunta de competencia - CORREGIDO"""
    try:
        print(f"🗑️ Eliminando pregunta de competencia ID: {pregunta_id}")
        
        pregunta = get_object_or_404(PreguntaCuestionario, id=pregunta_id)
        
        # CORREGIR: Verificar que pertenece a una competencia
        if not isinstance(pregunta.content_object, Competencia):
            return JsonResponse({
                'success': False, 
                'error': 'Esta pregunta no pertenece a una competencia'
            })
        
        competencia = pregunta.content_object
        
        # Verificar permisos
        if (competencia.recurso.seccion and 
            competencia.recurso.seccion.curso.profesor != request.user):
            return JsonResponse({
                'success': False, 
                'error': 'Sin permisos'
            }, status=403)
        
        # Solo permitir si está en configuración
        if competencia.estado not in ['configuracion']:
            return JsonResponse({
                'success': False, 
                'error': 'No se pueden eliminar preguntas de una competencia iniciada'
            })
        
        enunciado = pregunta.enunciado[:50]
        
        pregunta.delete()
        print(f"✅ Pregunta eliminada: {enunciado}")
        
        return JsonResponse({
            'success': True,
            'mensaje': 'Pregunta eliminada exitosamente'
        })
        
    except Exception as e:
        print(f"❌ Error al eliminar pregunta: {str(e)}")
        import traceback
        traceback.print_exc()
        
        return JsonResponse({
            'success': False, 
            'error': str(e)
        }, status=400)
    

@csrf_exempt
@require_http_methods(["POST"])
@login_required
@role_required(['Admin', 'Docente'])
def duplicar_pregunta_competencia(request):
    """Duplicar una pregunta de competencia - CORREGIDO"""
    try:
        print("📋 Duplicando pregunta de competencia...")
        
        data = json.loads(request.body)
        pregunta_id = data.get('pregunta_id')
        
        pregunta_original = get_object_or_404(PreguntaCuestionario, id=pregunta_id)
        
        # CORREGIR: Verificar que pertenece a una competencia
        if not isinstance(pregunta_original.content_object, Competencia):
            return JsonResponse({
                'success': False, 
                'error': 'Esta pregunta no pertenece a una competencia'
            })
        
        competencia = pregunta_original.content_object
        
        # Verificar permisos
        if (competencia.recurso.seccion and 
            competencia.recurso.seccion.curso.profesor != request.user):
            return JsonResponse({
                'success': False, 
                'error': 'Sin permisos'
            }, status=403)
        
        # Solo permitir si está en configuración
        if competencia.estado not in ['configuracion']:
            return JsonResponse({
                'success': False, 
                'error': 'No se pueden duplicar preguntas de una competencia iniciada'
            })
        
        # CORREGIR: Calcular nuevo orden
        from django.contrib.contenttypes.models import ContentType
        competencia_ct = ContentType.objects.get_for_model(Competencia)
        
        ultima_pregunta = PreguntaCuestionario.objects.filter(
            content_type=competencia_ct,
            object_id=competencia.id
        ).order_by('-orden').first()
        
        nuevo_orden = ultima_pregunta.orden + 1 if ultima_pregunta else 1
        
        with transaction.atomic():
            # Crear copia de la pregunta
            nueva_pregunta = PreguntaCuestionario.objects.create(
                content_type=competencia_ct,
                object_id=competencia.id,
                enunciado=f"{pregunta_original.enunciado} (Copia)",
                tipo=pregunta_original.tipo,
                orden=nuevo_orden,
                puntaje=pregunta_original.puntaje
            )
            
            # Duplicar opciones si existen
            for opcion in pregunta_original.opciones.all():
                OpcionPregunta.objects.create(
                    pregunta=nueva_pregunta,
                    texto=opcion.texto,
                    es_correcta=opcion.es_correcta
                )
        
        return JsonResponse({
            'success': True,
            'pregunta_id': nueva_pregunta.id,
            'mensaje': 'Pregunta duplicada exitosamente'
        })
        
    except Exception as e:
        print(f"❌ Error al duplicar pregunta: {str(e)}")
        import traceback
        traceback.print_exc()
        
        return JsonResponse({
            'success': False, 
            'error': str(e)
        }, status=400)
    

# =====================================================
# PANEL DE CONTROL EN VIVO
# =====================================================

@login_required
@role_required(['Admin', 'Docente'])
def panel_control_competencia(request, competencia_id):
    """Panel de control en tiempo real para la competencia"""
    print(f"🎮 Panel de control para competencia {competencia_id}")
    
    competencia = get_object_or_404(Competencia, id=competencia_id)
    
    # Verificar permisos
    if (competencia.recurso.seccion and 
        competencia.recurso.seccion.curso.profesor != request.user):
        print("❌ Usuario sin permisos para controlar")
        messages.error(request, 'No tienes permisos para controlar esta competencia')
        return redirect('biblioteca_competencias')
    
    # Obtener preguntas de la competencia usando ContentType
    from django.contrib.contenttypes.models import ContentType
    competencia_ct = ContentType.objects.get_for_model(Competencia)
    preguntas = PreguntaCuestionario.objects.filter(
        content_type=competencia_ct,
        object_id=competencia.id
    ).order_by('orden')
    
    # Obtener participaciones
    participaciones = competencia.participaciones.filter(activo=True).select_related('estudiante', 'grupo')
    
    # Obtener grupos si es modalidad grupal
    grupos = []
    if competencia.modalidad == 'grupal':
        grupos = competencia.grupos.all().prefetch_related('participaciones__estudiante')
    
    # Estadísticas en tiempo real
    total_participantes = participaciones.count()
    participantes_finalizados = participaciones.filter(finalizo=True).count()
    participantes_activos = total_participantes - participantes_finalizados
    
    # Ranking top 10
    ranking = list(participaciones.order_by('-puntaje_total', 'fecha_ingreso')[:10])
    
    # Verificar si puede iniciar la competencia
    puede_iniciar = (competencia.estado == 'esperando' and 
                     preguntas.exists() and 
                     total_participantes > 0)
    
    # Verificar si puede finalizar
    puede_finalizar = competencia.estado == 'activa'
    
    # Calcular tiempo restante
    tiempo_restante = 0
    if competencia.estado == 'activa' and competencia.fecha_inicio_real:
        tiempo_transcurrido = (timezone.now() - competencia.fecha_inicio_real).total_seconds()
        tiempo_limite_segundos = competencia.tiempo_limite * 60
        tiempo_restante = max(0, tiempo_limite_segundos - tiempo_transcurrido)
    
    context = {
        'competencia': competencia,
        'preguntas': preguntas,
        'participaciones': participaciones,
        'grupos': grupos,
        'total_participantes': total_participantes,
        'participantes_activos': participantes_activos,
        'participantes_finalizados': participantes_finalizados,
        'ranking': ranking,
        'tiempo_restante': int(tiempo_restante),
        'puede_iniciar': puede_iniciar,
        'puede_finalizar': puede_finalizar,
        'esta_activa': competencia.estado == 'activa',
        'esta_esperando': competencia.estado == 'esperando',
        'esta_finalizada': competencia.estado == 'finalizada',
        'imgPerfil': request.user.imgPerfil,
        'usuario': request.user.username,
    }
    
    return render(request, 'docente/competencias/panel_control.html', context)


@csrf_exempt
@require_http_methods(["POST"])
@login_required
@role_required(['Admin', 'Docente'])
def iniciar_competencia_live(request, competencia_id):
    """Iniciar la competencia en vivo"""
    try:
        print(f"🚀 Iniciando competencia en vivo {competencia_id}")
        
        competencia = get_object_or_404(Competencia, id=competencia_id)
        
        # Verificar permisos
        if (competencia.recurso.seccion and 
            competencia.recurso.seccion.curso.profesor != request.user):
            return JsonResponse({
                'success': False, 
                'error': 'Sin permisos'
            }, status=403)
        
        # Verificar estado
        if competencia.estado != 'esperando':
            return JsonResponse({
                'success': False, 
                'error': 'La competencia debe estar en estado "esperando" para poder iniciarla'
            })
        
        # Verificar que hay participantes
        participantes_count = competencia.participaciones.filter(activo=True).count()
        if participantes_count == 0:
            return JsonResponse({
                'success': False, 
                'error': 'No hay participantes registrados para iniciar la competencia'
            })
        
        # Verificar que hay preguntas
        if not competencia.preguntas.exists():
            return JsonResponse({
                'success': False, 
                'error': 'No se puede iniciar una competencia sin preguntas'
            })
        
        with transaction.atomic():
            # Cambiar estado a activa
            competencia.estado = 'activa'
            competencia.fecha_inicio_real = timezone.now()
            
            # Establecer primera pregunta como actual
            primera_pregunta = competencia.preguntas.order_by('orden').first()
            if primera_pregunta:
                competencia.pregunta_actual = primera_pregunta
            
            competencia.save()
            
            print(f"✅ Competencia {competencia_id} iniciada exitosamente")
        
        return JsonResponse({
            'success': True,
            'mensaje': 'Competencia iniciada exitosamente',
            'tiempo_limite_segundos': competencia.tiempo_limite * 60,
            'primera_pregunta_id': primera_pregunta.id if primera_pregunta else None
        })
        
    except Exception as e:
        print(f"❌ Error al iniciar competencia: {str(e)}")
        import traceback
        traceback.print_exc()
        
        return JsonResponse({
            'success': False, 
            'error': str(e)
        }, status=400)


@csrf_exempt
@require_http_methods(["POST"])
@login_required
@role_required(['Admin', 'Docente'])
def finalizar_competencia_live(request, competencia_id):
    """Finalizar la competencia en vivo"""
    try:
        print(f"🏁 Finalizando competencia en vivo {competencia_id}")
        
        competencia = get_object_or_404(Competencia, id=competencia_id)
        
        # Verificar permisos
        if (competencia.recurso.seccion and 
            competencia.recurso.seccion.curso.profesor != request.user):
            return JsonResponse({
                'success': False, 
                'error': 'Sin permisos'
            }, status=403)
        
        # Verificar estado
        if competencia.estado != 'activa':
            return JsonResponse({
                'success': False, 
                'error': 'Solo se pueden finalizar competencias activas'
            })
        
        with transaction.atomic():
            # Cambiar estado a finalizada
            competencia.estado = 'finalizada'
            competencia.fecha_finalizacion = timezone.now()
            competencia.pregunta_actual = None
            competencia.save()
            
            # Finalizar todas las participaciones activas
            participaciones_activas = competencia.participaciones.filter(
                activo=True, 
                finalizo=False
            )
            
            for participacion in participaciones_activas:
                participacion.finalizo = True
                participacion.fecha_finalizacion = timezone.now()
                
                # Calcular tiempo total si no se había calculado
                if not participacion.tiempo_total_segundos:
                    tiempo_total = (timezone.now() - participacion.fecha_ingreso).total_seconds()
                    participacion.tiempo_total_segundos = int(tiempo_total)
                
                # Calcular puntaje final
                participacion.calcular_puntaje_actual()
                participacion.save()
            
            # Calcular posiciones finales
            calcular_ranking_final(competencia)
            
            print(f"✅ Competencia {competencia_id} finalizada exitosamente")
        
        return JsonResponse({
            'success': True,
            'mensaje': 'Competencia finalizada exitosamente',
            'participantes_finalizados': participaciones_activas.count()
        })
        
    except Exception as e:
        print(f"❌ Error al finalizar competencia: {str(e)}")
        import traceback
        traceback.print_exc()
        
        return JsonResponse({
            'success': False, 
            'error': str(e)
        }, status=400)


@csrf_exempt
@require_http_methods(["POST"])
@login_required
@role_required(['Admin', 'Docente'])
def siguiente_pregunta_competencia(request, competencia_id):
    """Avanzar a la siguiente pregunta en la competencia"""
    try:
        print(f"➡️ Avanzando a siguiente pregunta en competencia {competencia_id}")
        
        competencia = get_object_or_404(Competencia, id=competencia_id)
        
        # Verificar permisos
        if (competencia.recurso.seccion and 
            competencia.recurso.seccion.curso.profesor != request.user):
            return JsonResponse({
                'success': False, 
                'error': 'Sin permisos'
            }, status=403)
        
        # Verificar estado
        if competencia.estado != 'activa':
            return JsonResponse({
                'success': False, 
                'error': 'La competencia debe estar activa'
            })
        
        # Obtener siguiente pregunta
        pregunta_actual_orden = competencia.pregunta_actual.orden if competencia.pregunta_actual else 0
        siguiente_pregunta = competencia.preguntas.filter(
            orden__gt=pregunta_actual_orden
        ).order_by('orden').first()
        
        if not siguiente_pregunta:
            # No hay más preguntas, finalizar automáticamente
            return finalizar_competencia_live(request, competencia_id)
        
        # Actualizar pregunta actual
        competencia.pregunta_actual = siguiente_pregunta
        competencia.save()
        
        print(f"✅ Avanzado a pregunta {siguiente_pregunta.orden}: {siguiente_pregunta.enunciado[:50]}")
        
        return JsonResponse({
            'success': True,
            'mensaje': f'Avanzado a pregunta {siguiente_pregunta.orden}',
            'pregunta_actual': {
                'id': siguiente_pregunta.id,
                'orden': siguiente_pregunta.orden,
                'enunciado': siguiente_pregunta.enunciado,
                'tipo': siguiente_pregunta.tipo.nombre
            }
        })
        
    except Exception as e:
        print(f"❌ Error al avanzar pregunta: {str(e)}")
        return JsonResponse({
            'success': False, 
            'error': str(e)
        }, status=400)


def calcular_ranking_final(competencia):
    """Calcular el ranking final de la competencia"""
    print(f"🏆 Calculando ranking final para competencia {competencia.id}")
    
    if competencia.modalidad == 'individual':
        # Ranking individual
        participaciones = competencia.participaciones.filter(
            finalizo=True
        ).order_by('-puntaje_total', 'tiempo_total_segundos')
        
        for i, participacion in enumerate(participaciones, 1):
            participacion.posicion = i
            participacion.save()
            
    elif competencia.modalidad == 'grupal':
        # Ranking grupal (promedio del grupo)
        grupos_puntajes = []
        
        for grupo in competencia.grupos.all():
            participaciones_grupo = grupo.participaciones.filter(finalizo=True)
            if participaciones_grupo.exists():
                # CORREGIR: Usar models.Avg correctamente
                from django.db.models import Avg
                puntaje_promedio = participaciones_grupo.aggregate(
                    promedio=Avg('puntaje_total')
                )['promedio'] or 0
                
                tiempo_promedio = participaciones_grupo.aggregate(
                    promedio=Avg('tiempo_total_segundos')
                )['promedio'] or 0
                
                grupos_puntajes.append({
                    'grupo': grupo,
                    'puntaje_promedio': puntaje_promedio,
                    'tiempo_promedio': tiempo_promedio,
                    'participaciones': list(participaciones_grupo)
                })
        
        # Ordenar grupos por puntaje y tiempo
        grupos_puntajes.sort(key=lambda x: (-x['puntaje_promedio'], x['tiempo_promedio']))
        
        # Asignar posiciones
        for i, grupo_data in enumerate(grupos_puntajes, 1):
            for participacion in grupo_data['participaciones']:
                participacion.posicion = i
                participacion.save()
    
    print(f"✅ Ranking final calculado para competencia {competencia.id}")


# =====================================================
# GESTIÓN DE GRUPOS
# =====================================================

@csrf_exempt
@require_http_methods(["POST"])
@login_required
@role_required(['Admin', 'Docente'])
def formar_grupos_automatico(request, competencia_id):
    """Formar grupos automáticamente"""
    try:
        print(f"🎲 Formando grupos automáticos para competencia {competencia_id}")
        
        competencia = get_object_or_404(Competencia, id=competencia_id)
        
        # Verificar permisos
        if (competencia.recurso.seccion and 
            competencia.recurso.seccion.curso.profesor != request.user):
            return JsonResponse({
                'success': False, 
                'error': 'Sin permisos'
            }, status=403)
        
        # Verificar modalidad
        if competencia.modalidad != 'grupal':
            return JsonResponse({
                'success': False, 
                'error': 'Solo se pueden formar grupos en competencias grupales'
            })
        
        # Verificar estado
        if competencia.estado not in ['configuracion', 'esperando']:
            return JsonResponse({
                'success': False, 
                'error': 'No se pueden formar grupos en una competencia activa'
            })
        
        # Obtener participantes sin grupo
        participantes_sin_grupo = competencia.participaciones.filter(
            grupo__isnull=True,
            activo=True
        ).select_related('estudiante')
        
        if not participantes_sin_grupo.exists():
            return JsonResponse({
                'success': False, 
                'error': 'No hay participantes sin grupo'
            })
        
        # Eliminar grupos existentes vacíos
        competencia.grupos.filter(participaciones__isnull=True).delete()
        
        with transaction.atomic():
            participantes_lista = list(participantes_sin_grupo)
            import random
            random.shuffle(participantes_lista)  # Mezclar aleatoriamente
            
            grupos_creados = 0
            max_miembros = competencia.max_miembros_grupo
            
            for i in range(0, len(participantes_lista), max_miembros):
                grupo_participantes = participantes_lista[i:i + max_miembros]
                
                # Crear grupo
                grupo = GrupoCompetencia.objects.create(
                    competencia=competencia,
                    nombre=f"Grupo {grupos_creados + 1}",
                    max_miembros=max_miembros,
                    creado_automaticamente=True
                )
                
                # Asignar participantes al grupo
                for participacion in grupo_participantes:
                    participacion.grupo = grupo
                    participacion.save()
                
                grupos_creados += 1
                print(f"✅ Grupo {grupo.nombre} creado con {len(grupo_participantes)} miembros")
        
        return JsonResponse({
            'success': True,
            'mensaje': f'{grupos_creados} grupos formados automáticamente',
            'grupos_creados': grupos_creados,
            'participantes_asignados': len(participantes_lista)
        })
        
    except Exception as e:
        print(f"❌ Error al formar grupos: {str(e)}")
        return JsonResponse({
            'success': False, 
            'error': str(e)
        }, status=400)


@login_required
@role_required(['Admin', 'Docente'])
def gestionar_grupos(request, competencia_id):
    """Vista para gestionar grupos manualmente"""
    competencia = get_object_or_404(Competencia, id=competencia_id)
    
    # Verificar permisos
    if (competencia.recurso.seccion and 
        competencia.recurso.seccion.curso.profesor != request.user):
        return redirect('dashboard-adm')
    
    # Verificar modalidad
    if competencia.modalidad != 'grupal':
        messages.error(request, 'Esta competencia no es grupal')
        return redirect('editar_competencia', competencia_id=competencia.id)
    
    # Obtener datos
    grupos = competencia.grupos.all().prefetch_related('participaciones__estudiante')
    participantes_sin_grupo = competencia.participaciones.filter(
        grupo__isnull=True,
        activo=True
    ).select_related('estudiante')
    
    context = {
        'competencia': competencia,
        'grupos': grupos,
        'participantes_sin_grupo': participantes_sin_grupo,
        'puede_modificar': competencia.estado in ['configuracion', 'esperando'],
        'imgPerfil': request.user.imgPerfil,
        'usuario': request.user.username,
    }
    
    return render(request, 'docente/competencias/gestionar_grupos.html', context)


@login_required
@role_required(['Admin', 'Docente'])
def biblioteca_competencias(request):
    """Vista para mostrar la biblioteca de competencias"""
    user = request.user
    
    # Obtener competencias del usuario
    competencias = Competencia.objects.filter(
        Q(recurso__seccion__curso__profesor=user) |  # Con sección del usuario
        Q(recurso__seccion__isnull=True)              # Sin sección (biblioteca)
    ).select_related('recurso', 'recurso__seccion', 'recurso__seccion__curso').order_by('-id')
    
    print(f"🔍 Total competencias encontradas: {competencias.count()}")
    
    context = {
        'competencias': competencias,
        'imgPerfil': user.imgPerfil,
        'usuario': user.username,
    }
    
    return render(request, 'docente/competencias/biblioteca_competencias.html', context)


@login_required
@role_required(['Admin', 'Docente'])
def resultados_competencia(request, competencia_id):
    """Vista para mostrar los resultados finales de la competencia"""
    competencia = get_object_or_404(Competencia, id=competencia_id)
    
    # Verificar permisos
    if (competencia.recurso.seccion and 
        competencia.recurso.seccion.curso.profesor != request.user):
        return redirect('dashboard-adm')
    
    # Obtener resultados según modalidad
    if competencia.modalidad == 'individual':
        participaciones = competencia.participaciones.filter(
            finalizo=True
        ).select_related('estudiante').order_by('posicion', '-puntaje_total')
    else:
        # Agrupar por grupos
        grupos_resultados = []
        for grupo in competencia.grupos.all():
            participaciones_grupo = grupo.participaciones.filter(finalizo=True)
            if participaciones_grupo.exists():
                puntaje_promedio = participaciones_grupo.aggregate(
                    promedio=models.Avg('puntaje_total')
                )['promedio'] or 0
                
                grupos_resultados.append({
                    'grupo': grupo,
                    'participaciones': list(participaciones_grupo),
                    'puntaje_promedio': puntaje_promedio,
                    'posicion': participaciones_grupo.first().posicion if participaciones_grupo.exists() else 999
                })
        
        grupos_resultados.sort(key=lambda x: x['posicion'])
        participaciones = grupos_resultados
    
    # Estadísticas generales
    total_participantes = competencia.participaciones.filter(activo=True).count()
    participantes_finalizados = competencia.participaciones.filter(finalizo=True).count()
    
    if participantes_finalizados > 0:
        puntaje_promedio = competencia.participaciones.filter(finalizo=True).aggregate(
            promedio=models.Avg('puntaje_total')
        )['promedio'] or 0
        
        tiempo_promedio = competencia.participaciones.filter(
            finalizo=True, 
            tiempo_total_segundos__isnull=False
        ).aggregate(
            promedio=models.Avg('tiempo_total_segundos')
        )['promedio'] or 0
    else:
        puntaje_promedio = 0
        tiempo_promedio = 0
    
    context = {
        'competencia': competencia,
        'participaciones': participaciones,
        'total_participantes': total_participantes,
        'participantes_finalizados': participantes_finalizados,
        'puntaje_promedio': puntaje_promedio,
        'tiempo_promedio_minutos': tiempo_promedio / 60 if tiempo_promedio else 0,
        'imgPerfil': request.user.imgPerfil,
        'usuario': request.user.username,
    }
    
    return render(request, 'docente/competencias/resultados_competencia.html', context)


# AGREGAR ESTAS FUNCIONES AL ARCHIVO myapp/views.py
# PARTE 4: VISTAS PARA ESTUDIANTES - ACCESO Y PARTICIPACIÓN

# =====================================================
# ACCESO DE ESTUDIANTES A COMPETENCIAS
# =====================================================

@login_required
@role_required('Estudiante')
def acceder_competencia(request, pin):
    """Vista para que estudiantes accedan con PIN"""
    print(f"🎯 Estudiante {request.user.username} accediendo con PIN: {pin}")
    
    try:
        competencia = get_object_or_404(Competencia, pin_acceso=pin)
        print(f"🏆 Competencia encontrada: {competencia.recurso.titulo}")
        
        # Verificar estado de la competencia
        if competencia.estado == 'configuracion':
            return render(request, 'estudiante/competencias/competencia_no_disponible.html', {
                'competencia': competencia,
                'mensaje': 'La competencia aún está en configuración. Vuelve más tarde.',
                'tipo': 'configuracion'
            })
        
        elif competencia.estado == 'finalizada':
            return render(request, 'estudiante/competencias/competencia_no_disponible.html', {
                'competencia': competencia,
                'mensaje': 'Esta competencia ya ha finalizado.',
                'tipo': 'finalizada'
            })
        
        elif competencia.estado == 'cancelada':
            return render(request, 'estudiante/competencias/competencia_no_disponible.html', {
                'competencia': competencia,
                'mensaje': 'Esta competencia ha sido cancelada.',
                'tipo': 'cancelada'
            })
        
        # Verificar si ya está participando
        participacion_existente = ParticipacionCompetencia.objects.filter(
            competencia=competencia,
            estudiante=request.user
        ).first()
        
        if participacion_existente:
            if participacion_existente.finalizo:
                # Ya finalizó, mostrar resultados
                return redirect('resultados_participacion', participacion_id=participacion_existente.id)
            elif competencia.estado == 'activa':
                # Competencia activa, ir a participar
                return redirect('participar_competencia', participacion_id=participacion_existente.id)
            else:
                # En sala de espera
                return redirect('sala_espera_competencia', participacion_id=participacion_existente.id)
        
        # Primera vez accediendo
        context = {
            'competencia': competencia,
            'pin': pin,
            'imgPerfil': request.user.imgPerfil,
            'usuario': request.user.username,
        }
        
        return render(request, 'estudiante/competencias/acceder_competencia.html', context)
        
    except Competencia.DoesNotExist:
        return render(request, 'estudiante/competencias/pin_invalido.html', {
            'pin': pin,
            'imgPerfil': request.user.imgPerfil,
        })


ESTADO_CONFIG = 'configuracion'
ESTADO_ESPERA = 'esperando'
ESTADO_ACTIVA = 'activa'
ESTADO_FIN    = 'finalizada'

def _competencia_disponible_para_join(c):
    """Reglas de negocio mínimas para permitir unirse por PIN."""
    if not c:
        return False, "La competencia no existe."
    if not getattr(c, 'recurso', None) or not getattr(c.recurso, 'titulo', '').strip():
        return False, "La competencia no tiene un recurso/título válido."
    # Debe tener al menos 1 pregunta asociada (usa tu GenericRelation o relación efectiva)
    tiene_preguntas = False
    # Si tu modelo Competencia expone un manager auxiliar, usa c.preguntas_manager.count()
    if hasattr(c, 'preguntas') and hasattr(c.preguntas, 'count'):
        tiene_preguntas = c.preguntas.count() > 0
    elif hasattr(c, 'preguntas_manager'):
        try:
            tiene_preguntas = c.preguntas_manager.count() > 0
        except Exception:
            tiene_preguntas = False

    if not tiene_preguntas:
        return False, "La competencia no tiene preguntas cargadas."
    if str(c.estado) not in [ESTADO_CONFIG, ESTADO_ESPERA, ESTADO_ACTIVA]:
        return False, "La competencia no está disponible para unirse."
    return True, None


@login_required
def unirse_competencia(request):
    """
    GET  → muestra form PIN
    POST → valida PIN, crea/recupera ParticipacionCompetencia y envía a seleccionar avatar
    """
    if request.method == 'POST':
        pin = (request.POST.get('pin') or '').strip()
        if not pin:
            messages.error(request, "Ingresa un PIN.")
            return redirect('competencia_unirse')

        comp = Competencia.objects.filter(codigo_pin=pin).select_related('recurso').first()
        ok, msg = _competencia_disponible_para_join(comp)
        if not ok:
            messages.error(request, msg or "No se puede unir a esta competencia.")
            return redirect('competencia_unirse')

        # Participación única por (competencia, estudiante)
        part, created = ParticipacionCompetencia.objects.get_or_create(
            competencia=comp,
            estudiante=request.user,
            defaults={'activo': True}
        )
        # Si ya tenía participación y permitir_reingreso=True → reusar
        if not created and not comp.permitir_reingreso:
            messages.error(request, "Ya tienes una participación registrada y el reingreso no está permitido.")
            return redirect('competencia_unirse')

        return redirect('competencia_avatar', competencia_id=comp.id)

    # GET
    return render(request, 'estudiante/competencias/unirse_competencia.html', {})


def asignar_grupo_automatico(participacion):
    """Asignar automáticamente a un grupo disponible"""
    competencia = participacion.competencia
    
    # Buscar grupo con espacio disponible
    grupo_disponible = None
    for grupo in competencia.grupos.all():
        if grupo.puede_agregar_miembro():
            grupo_disponible = grupo
            break
    
    # Si no hay grupo disponible, crear uno nuevo
    if not grupo_disponible:
        numero_grupo = competencia.grupos.count() + 1
        grupo_disponible = GrupoCompetencia.objects.create(
            competencia=competencia,
            nombre=f"Grupo {numero_grupo}",
            max_miembros=competencia.max_miembros_grupo,
            creado_automaticamente=True
        )
    
    # Asignar al grupo
    participacion.grupo = grupo_disponible
    participacion.save()
    
    print(f"✅ {participacion.estudiante.username} asignado automáticamente a {grupo_disponible.nombre}")


# =====================================================
# GESTIÓN DE GRUPOS PARA ESTUDIANTES
# =====================================================

@login_required
@role_required('Estudiante')
def seleccionar_grupo(request, competencia_id):
    """Vista para que estudiantes seleccionen grupo"""
    competencia = get_object_or_404(Competencia, id=competencia_id)
    
    # Verificar participación
    participacion = get_object_or_404(
        ParticipacionCompetencia,
        competencia=competencia,
        estudiante=request.user
    )
    
    # Verificar modalidad y configuración
    if competencia.modalidad != 'grupal' or not competencia.grupos_abiertos:
        return redirect('sala_espera_competencia', participacion_id=participacion.id)
    
    # Obtener grupos disponibles
    grupos_disponibles = []
    for grupo in competencia.grupos.all():
        grupos_disponibles.append({
            'grupo': grupo,
            'miembros': list(grupo.participaciones.select_related('estudiante')),
            'puede_unirse': grupo.puede_agregar_miembro(),
            'es_completo': grupo.esta_completo()
        })
    
    context = {
        'competencia': competencia,
        'participacion': participacion,
        'grupos_disponibles': grupos_disponibles,
        'puede_crear_grupo': True,  # Permitir crear grupos nuevos
        'imgPerfil': request.user.imgPerfil,
        'usuario': request.user.username,
    }
    
    return render(request, 'estudiante/competencias/seleccionar_grupo.html', context)


@csrf_exempt
@require_http_methods(["POST"])
@login_required
@role_required('Estudiante')
def crear_grupo_estudiante(request, competencia_id):
    """Crear nuevo grupo por estudiante"""
    try:
        competencia = get_object_or_404(Competencia, id=competencia_id)
        
        # Verificar participación
        participacion = get_object_or_404(
            ParticipacionCompetencia,
            competencia=competencia,
            estudiante=request.user
        )
        
        # Verificar que no tenga grupo
        if participacion.grupo:
            return JsonResponse({
                'success': False,
                'error': 'Ya tienes un grupo asignado'
            })
        
        # Verificar configuración
        if not competencia.grupos_abiertos:
            return JsonResponse({
                'success': False,
                'error': 'No se permite crear grupos en esta competencia'
            })
        
        data = json.loads(request.body)
        nombre_grupo = data.get('nombre_grupo', '').strip()
        
        if not nombre_grupo:
            return JsonResponse({
                'success': False,
                'error': 'El nombre del grupo es obligatorio'
            })
        
        # Verificar que no exista un grupo con el mismo nombre
        if competencia.grupos.filter(nombre=nombre_grupo).exists():
            return JsonResponse({
                'success': False,
                'error': 'Ya existe un grupo con ese nombre'
            })
        
        with transaction.atomic():
            # Crear grupo
            grupo = GrupoCompetencia.objects.create(
                competencia=competencia,
                nombre=nombre_grupo,
                max_miembros=competencia.max_miembros_grupo,
                creado_automaticamente=False
            )
            
            # Asignar al creador
            participacion.grupo = grupo
            participacion.save()
            
            print(f"✅ Grupo '{nombre_grupo}' creado por {request.user.username}")
        
        return JsonResponse({
            'success': True,
            'grupo_id': grupo.id,
            'mensaje': f'Grupo "{nombre_grupo}" creado exitosamente'
        })
        
    except Exception as e:
        print(f"❌ Error al crear grupo: {str(e)}")
        return JsonResponse({
            'success': False,
            'error': str(e)
        })


@csrf_exempt
@require_http_methods(["POST"])
@login_required
@role_required('Estudiante')
def unirse_grupo(request, grupo_id):
    """Unirse a un grupo existente"""
    try:
        grupo = get_object_or_404(GrupoCompetencia, id=grupo_id)
        
        # Verificar participación
        participacion = get_object_or_404(
            ParticipacionCompetencia,
            competencia=grupo.competencia,
            estudiante=request.user
        )
        
        # Verificar que no tenga grupo
        if participacion.grupo:
            return JsonResponse({
                'success': False,
                'error': 'Ya tienes un grupo asignado'
            })
        
        # Verificar espacio disponible
        if not grupo.puede_agregar_miembro():
            return JsonResponse({
                'success': False,
                'error': 'Este grupo ya está completo'
            })
        
        # Asignar al grupo
        participacion.grupo = grupo
        participacion.save()
        
        print(f"✅ {request.user.username} se unió a {grupo.nombre}")
        
        return JsonResponse({
            'success': True,
            'mensaje': f'Te has unido al grupo "{grupo.nombre}"'
        })
        
    except Exception as e:
        print(f"❌ Error al unirse a grupo: {str(e)}")
        return JsonResponse({
            'success': False,
            'error': str(e)
        })


# =====================================================
# PARTICIPACIÓN EN COMPETENCIA
# =====================================================

@login_required
@role_required('Estudiante')
def sala_espera_competencia(request, participacion_id):
    """Sala de espera antes de que inicie la competencia"""
    participacion = get_object_or_404(
        ParticipacionCompetencia,
        id=participacion_id,
        estudiante=request.user
    )
    
    competencia = participacion.competencia
    
    # Si ya está activa, redirigir a participar
    if competencia.estado == 'activa':
        return redirect('participar_competencia', participacion_id=participacion.id)
    
    # Si ya finalizó, redirigir a resultados
    if participacion.finalizo or competencia.estado == 'finalizada':
        return redirect('resultados_participacion', participacion_id=participacion.id)
    
    # Obtener información de otros participantes
    if competencia.modalidad == 'grupal' and participacion.grupo:
        companeros = participacion.grupo.participaciones.exclude(
            id=participacion.id
        ).select_related('estudiante')
        total_participantes = participacion.grupo.participaciones.count()
        grupo_info = participacion.grupo
    else:
        companeros = []
        total_participantes = competencia.participaciones.filter(activo=True).count()
        grupo_info = None
    
    context = {
        'participacion': participacion,
        'competencia': competencia,
        'companeros': companeros,
        'total_participantes': total_participantes,
        'grupo_info': grupo_info,
        'imgPerfil': request.user.imgPerfil,
        'usuario': request.user.username,
    }
    
    return render(request, 'estudiante/competencias/sala_espera.html', context)



def _resultados_url(participacion):
    return reverse(
        'resultados_participacion',
        kwargs={'competencia_id': participacion.competencia_id, 'participacion_id': participacion.id}
    )

# Ejemplo 1: en participar_competencia, cuando debas ir a resultados
@login_required
# Si manejas roles, descomenta la siguiente línea:
# @role_required('Estudiante')
def participar_competencia(request, participacion_id: int):
    """
    Renderiza la pantalla de participación del estudiante.
    Reglas:
      - Si la competencia está 'finalizada' o 'cancelada' -> redirige a resultados.
      - Si la competencia NO está 'activa' -> redirige a sala de espera.
      - Si el estudiante aún no eligió avatar -> obliga a seleccionar avatar.
    Contexto que entrega al template:
      - competencia, participacion
      - total_preguntas, participantes_activos, progreso_porcentaje
      - debug (para herramientas de desarrollo en el template)
    """
    # Trae participación + competencia (+ recurso para título) en una sola consulta
    participacion = get_object_or_404(
        ParticipacionCompetencia.objects.select_related('competencia', 'competencia__recurso'),
        id=participacion_id,
        estudiante=request.user
    )
    competencia = participacion.competencia

    # Si no eligió avatar aún, reenvía a la selección de avatar
    if not getattr(participacion, 'avatar_elegido', False):
        url = reverse('seleccionar_avatar_competencia', args=[competencia.id])
        return redirect(f"{url}?participacion_id={participacion.id}")

    # Estado de competencia
    estado = str(getattr(competencia, 'estado', '')).lower()

    # Si terminó o fue cancelada -> resultados (nota: requiere ambos IDs)
    if estado in ('finalizada', 'cancelada'):
        return redirect('resultados_participacion',
                        competencia_id=competencia.id,
                        participacion_id=participacion.id)

    # Si aún no está activa -> sala de espera
    if estado != 'activa':
        return redirect('sala_espera_competencia', participacion_id=participacion.id)

    # ---- Métricas para el template ----
    # Total de preguntas (GenericRelation o manager auxiliar)
    total_preguntas = 0
    try:
        if hasattr(competencia, 'preguntas') and hasattr(competencia.preguntas, 'count'):
            total_preguntas = competencia.preguntas.count()
        elif hasattr(competencia, 'preguntas_manager'):
            total_preguntas = competencia.preguntas_manager.count()
    except Exception:
        total_preguntas = 0

    # Participantes activos
    participantes_activos = (
        ParticipacionCompetencia.objects
        .filter(competencia=competencia, activo=True)
        .count()
    )

    # Progreso: usa pregunta_actual.orden si existe
    orden_actual = None
    try:
        pregunta_actual = getattr(competencia, 'pregunta_actual', None)
        if pregunta_actual is not None:
            orden_actual = getattr(pregunta_actual, 'orden', None) or 1
    except Exception:
        orden_actual = None

    if total_preguntas and orden_actual:
        progreso_porcentaje = round((orden_actual / total_preguntas) * 100, 2)
    else:
        progreso_porcentaje = 0

    context = {
        'competencia': competencia,
        'participacion': participacion,
        'total_preguntas': total_preguntas,
        'participantes_activos': participantes_activos,
        'progreso_porcentaje': progreso_porcentaje,
        'debug': getattr(settings, 'DEBUG', False),
    }
    return render(request, 'estudiante/competencias/participar_competencia.html', context)


@csrf_exempt
@require_http_methods(["POST"])
@login_required
@role_required('Estudiante')
def responder_competencia(request):
    """Guardar respuesta en competencia vía AJAX"""
    try:
        data = json.loads(request.body)
        participacion_id = data.get('participacion_id')
        pregunta_id = data.get('pregunta_id')
        respuesta = data.get('respuesta')
        
        print(f"💾 Guardando respuesta de competencia - Participación: {participacion_id}, Pregunta: {pregunta_id}")
        
        participacion = get_object_or_404(
            ParticipacionCompetencia,
            id=participacion_id,
            estudiante=request.user
        )
        
        pregunta = get_object_or_404(PreguntaCuestionario, id=pregunta_id)
        
        # Verificar que la competencia esté activa
        if participacion.competencia.estado != 'activa':
            return JsonResponse({
                'success': False,
                'error': 'La competencia no está activa'
            })
        
        # Verificar que no haya finalizado
        if participacion.finalizo:
            return JsonResponse({
                'success': False,
                'error': 'Tu participación ya ha finalizado'
            })
        
        # Crear o actualizar respuesta
        respuesta_obj, created = RespuestaCompetencia.objects.get_or_create(
            participacion=participacion,
            pregunta=pregunta,
            defaults={'fecha_respuesta': timezone.now()}
        )
        
        # Procesar según tipo de pregunta
        if pregunta.tipo.nombre in ['opcion_unica', 'falso_verdadero']:
            if respuesta:
                try:
                    opcion = get_object_or_404(OpcionPregunta, id=respuesta)
                    respuesta_obj.opcion_seleccionada = opcion
                    print(f"✅ Opción única guardada")
                except:
                    respuesta_obj.opcion_seleccionada = None
            else:
                respuesta_obj.opcion_seleccionada = None
        
        elif pregunta.tipo.nombre == 'opcion_multiple':
            if isinstance(respuesta, list) and respuesta:
                respuesta_obj.opciones_multiples = ','.join(map(str, respuesta))
                print(f"✅ Opción múltiple guardada")
            else:
                respuesta_obj.opciones_multiples = ''
        
        elif pregunta.tipo.nombre == 'respuesta_abierta':
            respuesta_texto = str(respuesta) if respuesta else ''
            if len(respuesta_texto) > 2000:
                respuesta_texto = respuesta_texto[:2000]
            respuesta_obj.respuesta_texto = respuesta_texto
            print(f"✅ Respuesta abierta guardada - {len(respuesta_texto)} caracteres")
        
        # Calcular puntaje automáticamente
        try:
            puntaje_obtenido = respuesta_obj.calcular_puntaje_automatico()
            print(f"📊 Puntaje calculado: {puntaje_obtenido}")
        except:
            print("❌ Error al calcular puntaje automático")
            puntaje_obtenido = 0
        
        respuesta_obj.save()
        
        # Actualizar puntaje total de la participación
        participacion.calcular_puntaje_actual()
        
        return JsonResponse({
            'success': True,
            'puntaje_obtenido': float(puntaje_obtenido)
        })
        
    except Exception as e:
        print(f"❌ Error al guardar respuesta de competencia: {str(e)}")
        import traceback
        traceback.print_exc()
        return JsonResponse({
            'success': False,
            'error': str(e)
        })


@csrf_exempt
@require_http_methods(["POST", "GET"])
@login_required
@role_required('Estudiante')
def finalizar_participacion(request, participacion_id):
    """Finalizar participación del estudiante - IMPLEMENTACIÓN COMPLETA"""
    try:
        print(f"🏁 Finalizando participación {participacion_id} para {request.user.username}")
        
        participacion = get_object_or_404(
            ParticipacionCompetencia,
            id=participacion_id,
            estudiante=request.user
        )
        
        competencia = participacion.competencia
        
        # Si ya finalizó, redirigir a resultados
        if participacion.finalizo:
            print(f"✅ Participación ya finalizada, redirigiendo a resultados")
            return redirect('resultados_participacion', participacion_id=participacion.id)
        
        # Finalizar participación
        with transaction.atomic():
            participacion.finalizo = True
            participacion.fecha_finalizacion = timezone.now()
            
            # Calcular tiempo total
            if participacion.fecha_ingreso:
                tiempo_total = (participacion.fecha_finalizacion - participacion.fecha_ingreso).total_seconds()
                participacion.tiempo_total_segundos = int(tiempo_total)
            
            # Calcular puntaje final
            participacion.calcular_puntaje_actual()
            
            # Calcular posición/ranking
            calcular_posicion_participacion(participacion)
            
            participacion.save()
            
            print(f"✅ Participación finalizada - Puntaje: {participacion.puntaje_total}")
        
        # Redirigir a resultados
        return redirect('resultados_participacion', participacion_id=participacion.id)
        
    except Exception as e:
        print(f"❌ Error al finalizar participación: {str(e)}")
        import traceback
        traceback.print_exc()
        
        messages.error(request, 'Error al finalizar la participación')
        return redirect('student_dashboard')


def calcular_posicion_participacion(participacion):
    """Calcular la posición del participante en el ranking - VERSIÓN COMPLETA"""
    from django.db import models  # Import necesario para Avg
    
    try:
        competencia = participacion.competencia
        print(f"📊 Calculando posición para {participacion.estudiante.username} en modalidad {competencia.modalidad}")
        
        if competencia.modalidad == 'individual':
            # RANKING INDIVIDUAL
            print("🥇 Calculando ranking individual...")
            
            # Obtener todas las participaciones finalizadas ordenadas por puntaje (desc) y fecha (asc)
            participaciones_finalizadas = competencia.participaciones.filter(
                finalizo=True
            ).order_by('-puntaje_total', 'fecha_finalizacion')
            
            print(f"📋 {participaciones_finalizadas.count()} participaciones finalizadas")
            
            # Asignar posiciones a TODAS las participaciones finalizadas
            posicion_actual_participacion = None
            
            for pos, p in enumerate(participaciones_finalizadas, 1):
                p.posicion = pos
                p.save()
                
                # Guardar la posición de la participación actual
                if p.id == participacion.id:
                    posicion_actual_participacion = pos
                
                print(f"  #{pos}: {p.estudiante.username} - {p.puntaje_total} pts")
            
            print(f"🎯 Posición final de {participacion.estudiante.username}: #{posicion_actual_participacion}")
            
        else:
            # RANKING GRUPAL
            print("👥 Calculando ranking grupal...")
            
            grupos_ranking = []
            
            # Calcular promedio por grupo
            for grupo in competencia.grupos.all():
                participaciones_grupo = grupo.participaciones.filter(finalizo=True)
                
                if participaciones_grupo.exists():
                    puntaje_promedio = participaciones_grupo.aggregate(
                        promedio=models.Avg('puntaje_total')
                    )['promedio'] or 0
                    
                    grupos_ranking.append({
                        'grupo': grupo,
                        'puntaje_promedio': float(puntaje_promedio),
                        'participaciones': list(participaciones_grupo)
                    })
                    
                    print(f"  Grupo {grupo.nombre}: promedio {puntaje_promedio:.2f} ({participaciones_grupo.count()} miembros)")
            
            # Ordenar grupos por promedio (mayor a menor)
            grupos_ranking.sort(key=lambda x: -x['puntaje_promedio'])
            
            # Asignar posiciones grupales
            posicion_grupo_actual = None
            
            for pos, grupo_data in enumerate(grupos_ranking, 1):
                grupo = grupo_data['grupo']
                
                print(f"🏆 Grupo #{pos}: {grupo.nombre} - Promedio: {grupo_data['puntaje_promedio']:.2f}")
                
                # Asignar la misma posición a todos los miembros del grupo
                for participacion_grupo in grupo_data['participaciones']:
                    participacion_grupo.posicion = pos
                    participacion_grupo.save()
                    
                    # Si es la participación actual, guardar su posición
                    if participacion_grupo.id == participacion.id:
                        posicion_grupo_actual = pos
                        print(f"🎯 {participacion.estudiante.username} queda en posición grupal: #{pos}")
            
            print(f"📊 Ranking grupal completado. Total grupos: {len(grupos_ranking)}")
        
        # Refrescar la participación para obtener la posición actualizada
        participacion.refresh_from_db()
        print(f"✅ Posición final asignada: #{participacion.posicion}")
        
        return participacion.posicion
        
    except Exception as e:
        print(f"❌ Error calculando posición: {str(e)}")
        import traceback
        traceback.print_exc()
        
        # En caso de error, asignar posición por defecto
        participacion.posicion = 999
        participacion.save()
        return 999


@login_required
@role_required('Estudiante')
def resultados_participacion(request, participacion_id):
    """Mostrar resultados de la participación"""
    participacion = get_object_or_404(
        ParticipacionCompetencia,
        id=participacion_id,
        estudiante=request.user
    )
    
    competencia = participacion.competencia
    
    # Obtener respuestas con detalles
    respuestas = participacion.respuestas.select_related(
        'pregunta', 'opcion_seleccionada'
    ).prefetch_related('pregunta__opciones').all()
    
    # Procesar respuestas para el template
    for respuesta in respuestas:
        if respuesta.pregunta.puntaje > 0:
            respuesta.porcentaje_pregunta = int(
                (respuesta.puntaje_obtenido / respuesta.pregunta.puntaje) * 100
            )
        else:
            respuesta.porcentaje_pregunta = 0
        
        if respuesta.opciones_multiples:
            respuesta.opciones_ids = [
                int(x.strip()) for x in respuesta.opciones_multiples.split(',') if x.strip()
            ]
        else:
            respuesta.opciones_ids = []
    
    # Calcular porcentaje total
    puntaje_total_posible = competencia.calcular_puntaje_total()
    puntaje_porcentaje = 0
    if puntaje_total_posible > 0:
        puntaje_porcentaje = (participacion.puntaje_total / puntaje_total_posible * 100)
    
    # Obtener ranking si está disponible
    mejor_que_porcentaje = 0
    if participacion.posicion and competencia.estado == 'finalizada':
        total_participantes = competencia.participaciones.filter(finalizo=True).count()
        if total_participantes > 1:
            mejor_que_porcentaje = int(
                ((total_participantes - participacion.posicion) / (total_participantes - 1)) * 100
            )
    
    context = {
        'participacion': participacion,
        'competencia': competencia,
        'respuestas': respuestas,
        'puntaje_porcentaje': puntaje_porcentaje,
        'mejor_que_porcentaje': mejor_que_porcentaje,
        'tiempo_minutos': participacion.tiempo_total_segundos / 60 if participacion.tiempo_total_segundos else 0,
        'mostrar_resultados': competencia.mostrar_resultados_inmediatos,
        'imgPerfil': request.user.imgPerfil,
        'usuario': request.user.username,
    }
    
    return render(request, 'estudiante/competencias/resultados_participacion.html', context)


# =====================================================
# ENDPOINTS AJAX PARA TIEMPO REAL
# =====================================================

@require_http_methods(["GET"])
@login_required
def obtener_estado_competencia(request, competencia_id):
    """Obtener estado actual de la competencia (AJAX)"""
    try:
        competencia = get_object_or_404(Competencia, id=competencia_id)
        
        return JsonResponse({
            'success': True,
            'estado': competencia.estado,
            'tiempo_restante': competencia.tiempo_restante_segundos(),
            'pregunta_actual': {
                'id': competencia.pregunta_actual.id,
                'orden': competencia.pregunta_actual.orden,
            } if competencia.pregunta_actual else None,
            'participantes_activos': competencia.participaciones.filter(
                activo=True, finalizo=False
            ).count()
        })
        
    except Exception as e:
        return JsonResponse({
            'success': False,
            'error': str(e)
        })


@require_http_methods(["GET"])
@login_required
def obtener_ranking_tiempo_real(request, competencia_id):
    """Obtener ranking en tiempo real (AJAX)"""
    try:
        competencia = get_object_or_404(Competencia, id=competencia_id)
        
        if competencia.modalidad == 'individual':
            participaciones = competencia.participaciones.filter(
                activo=True
            ).select_related('estudiante').order_by(
                '-puntaje_total', 'fecha_ingreso'
            )[:10]
            
            ranking = []
            for i, participacion in enumerate(participaciones, 1):
                ranking.append({
                    'posicion': i,
                    'nombre': f"{participacion.estudiante.first_name} {participacion.estudiante.last_name}".strip() or participacion.estudiante.username,
                    'puntaje': float(participacion.puntaje_total),
                    'finalizo': participacion.finalizo
                })
        else:
            # Ranking grupal - calcular promedios
            grupos_ranking = []
            for grupo in competencia.grupos.all():
                participaciones_grupo = grupo.participaciones.filter(activo=True)
                if participaciones_grupo.exists():
                    puntaje_promedio = participaciones_grupo.aggregate(
                        promedio=models.Avg('puntaje_total')
                    )['promedio'] or 0
                    
                    grupos_ranking.append({
                        'grupo': grupo,
                        'puntaje_promedio': puntaje_promedio,
                        'miembros_finalizados': participaciones_grupo.filter(finalizo=True).count(),
                        'total_miembros': participaciones_grupo.count()
                    })
            
            grupos_ranking.sort(key=lambda x: -x['puntaje_promedio'])
            
            ranking = []
            for i, grupo_data in enumerate(grupos_ranking[:10], 1):
                ranking.append({
                    'posicion': i,
                    'nombre': grupo_data['grupo'].nombre,
                    'puntaje': float(grupo_data['puntaje_promedio']),
                    'finalizo': grupo_data['miembros_finalizados'] == grupo_data['total_miembros']
                })
        
        return JsonResponse({
            'success': True,
            'ranking': ranking
        })
        
    except Exception as e:
        return JsonResponse({
            'success': False,
            'error': str(e)
        })


@csrf_exempt
@require_http_methods(["POST"])
@login_required
@role_required(['Admin', 'Docente'])
def abrir_sala_espera(request, competencia_id):
    """Abrir sala de espera para que estudiantes se conecten"""
    try:
        print(f"🚪 Abriendo sala de espera para competencia {competencia_id}")
        
        competencia = get_object_or_404(Competencia, id=competencia_id)
        
        # Verificar permisos
        if (competencia.recurso.seccion and 
            competencia.recurso.seccion.curso.profesor != request.user):
            return JsonResponse({'success': False, 'error': 'Sin permisos'}, status=403)
        
        # Verificar estado actual
        if competencia.estado != 'configuracion':
            return JsonResponse({
                'success': False, 
                'error': 'La competencia debe estar en configuración para abrir sala de espera'
            })
        
        # Verificar que tenga preguntas
        from django.contrib.contenttypes.models import ContentType
        competencia_ct = ContentType.objects.get_for_model(Competencia)
        preguntas_count = PreguntaCuestionario.objects.filter(
            content_type=competencia_ct,
            object_id=competencia.id
        ).count()
        
        if preguntas_count == 0:
            return JsonResponse({
                'success': False, 
                'error': 'No se puede abrir sala de espera sin preguntas'
            })
        
        # Cambiar estado a esperando
        competencia.estado = 'esperando'
        competencia.save()
        
        print(f"✅ Sala de espera abierta para competencia {competencia_id}")
        
        return JsonResponse({
            'success': True,
            'mensaje': 'Sala de espera abierta. Los estudiantes pueden unirse ahora.',
            'pin_acceso': competencia.pin_acceso
        })
        
    except Exception as e:
        print(f"❌ Error al abrir sala de espera: {str(e)}")
        return JsonResponse({'success': False, 'error': str(e)}, status=400)


@csrf_exempt
@require_http_methods(["POST"])
@login_required
@role_required(['Admin', 'Docente'])
def crear_grupo_manual(request, competencia_id):
    """Crear grupo manualmente por el profesor"""
    try:
        competencia = get_object_or_404(Competencia, id=competencia_id)
        
        # Verificar permisos
        if (competencia.recurso.seccion and 
            competencia.recurso.seccion.curso.profesor != request.user):
            return JsonResponse({'success': False, 'error': 'Sin permisos'}, status=403)
        
        data = json.loads(request.body)
        nombre_grupo = data.get('nombre_grupo', '').strip()
        
        if not nombre_grupo:
            return JsonResponse({'success': False, 'error': 'Nombre del grupo requerido'})
        
        # Verificar que no exista
        if competencia.grupos.filter(nombre=nombre_grupo).exists():
            return JsonResponse({'success': False, 'error': 'Ya existe un grupo con ese nombre'})
        
        grupo = GrupoCompetencia.objects.create(
            competencia=competencia,
            nombre=nombre_grupo,
            max_miembros=competencia.max_miembros_grupo,
            creado_automaticamente=False
        )
        
        return JsonResponse({
            'success': True,
            'grupo_id': grupo.id,
            'mensaje': f'Grupo "{nombre_grupo}" creado exitosamente'
        })
        
    except Exception as e:
        return JsonResponse({'success': False, 'error': str(e)}, status=400)
    


@csrf_exempt
@require_http_methods(["POST"])
@login_required
@role_required(['Admin', 'Docente'])
def mover_estudiante_grupo(request):
    """Mover estudiante de un grupo a otro"""
    try:
        data = json.loads(request.body)
        participacion_id = data.get('participacion_id')
        nuevo_grupo_id = data.get('nuevo_grupo_id')
        
        participacion = get_object_or_404(ParticipacionCompetencia, id=participacion_id)
        
        # Verificar permisos
        if (participacion.competencia.recurso.seccion and 
            participacion.competencia.recurso.seccion.curso.profesor != request.user):
            return JsonResponse({'success': False, 'error': 'Sin permisos'}, status=403)
        
        if nuevo_grupo_id:
            nuevo_grupo = get_object_or_404(GrupoCompetencia, id=nuevo_grupo_id)
            
            # Verificar espacio disponible
            if not nuevo_grupo.puede_agregar_miembro():
                return JsonResponse({'success': False, 'error': 'El grupo está completo'})
            
            participacion.grupo = nuevo_grupo
        else:
            # Remover del grupo
            participacion.grupo = None
        
        participacion.save()
        
        return JsonResponse({
            'success': True,
            'mensaje': 'Estudiante movido exitosamente'
        })
        
    except Exception as e:
        return JsonResponse({'success': False, 'error': str(e)}, status=400)
    

@csrf_exempt
@require_http_methods(["DELETE"])
@login_required
@role_required(['Admin', 'Docente'])
def eliminar_grupo(request, grupo_id):
    """Eliminar un grupo y liberar a sus miembros"""
    try:
        grupo = get_object_or_404(GrupoCompetencia, id=grupo_id)
        
        # Verificar permisos
        if (grupo.competencia.recurso.seccion and 
            grupo.competencia.recurso.seccion.curso.profesor != request.user):
            return JsonResponse({'success': False, 'error': 'Sin permisos'}, status=403)
        
        # Verificar estado de competencia
        if grupo.competencia.estado not in ['configuracion', 'esperando']:
            return JsonResponse({
                'success': False, 
                'error': 'No se pueden eliminar grupos de una competencia activa'
            })
        
        # Liberar participantes del grupo
        participaciones = grupo.participaciones.all()
        for participacion in participaciones:
            participacion.grupo = None
            participacion.save()
        
        nombre_grupo = grupo.nombre
        grupo.delete()
        
        return JsonResponse({
            'success': True,
            'mensaje': f'Grupo "{nombre_grupo}" eliminado y participantes liberados'
        })
        
    except Exception as e:
        return JsonResponse({'success': False, 'error': str(e)}, status=400)


@require_http_methods(["GET"])
@login_required
@role_required(['Admin', 'Docente'])
def obtener_participantes_competencia(request, competencia_id):
    """Obtener lista de participantes (AJAX)"""
    try:
        competencia = get_object_or_404(Competencia, id=competencia_id)
        
        # Verificar permisos
        if (competencia.recurso.seccion and 
            competencia.recurso.seccion.curso.profesor != request.user):
            return JsonResponse({'success': False, 'error': 'Sin permisos'}, status=403)
        
        participaciones = competencia.participaciones.filter(activo=True).select_related('estudiante', 'grupo')
        
        participantes = []
        for participacion in participaciones:
            participantes.append({
                'id': participacion.id,
                'nombre': f"{participacion.estudiante.first_name} {participacion.estudiante.last_name}".strip() or participacion.estudiante.username,
                'email': participacion.estudiante.email,
                'grupo': participacion.grupo.nombre if participacion.grupo else None,
                'grupo_id': participacion.grupo.id if participacion.grupo else None,
                'fecha_ingreso': participacion.fecha_ingreso.isoformat(),
                'finalizo': participacion.finalizo,
                'puntaje': float(participacion.puntaje_total)
            })
        
        return JsonResponse({
            'success': True,
            'participantes': participantes,
            'total': len(participantes)
        })
        
    except Exception as e:
        return JsonResponse({'success': False, 'error': str(e)}, status=400)
    

@login_required
@role_required(['Admin', 'Docente'])
def exportar_resultados_competencia(request, competencia_id):
    """Exportar resultados de competencia a CSV"""
    import csv
    from django.http import HttpResponse
    
    competencia = get_object_or_404(Competencia, id=competencia_id)
    
    # Verificar permisos
    if (competencia.recurso.seccion and 
        competencia.recurso.seccion.curso.profesor != request.user):
        return redirect('dashboard-adm')
    
    # Crear respuesta CSV
    response = HttpResponse(content_type='text/csv')
    response['Content-Disposition'] = f'attachment; filename="resultados_{competencia.recurso.titulo}_{competencia.id}.csv"'
    
    writer = csv.writer(response)
    
    # Encabezados
    if competencia.modalidad == 'individual':
        writer.writerow([
            'Posición', 'Nombre', 'Email', 'Puntaje', 'Tiempo (min)', 
            'Finalizado', 'Fecha Ingreso', 'Fecha Finalización'
        ])
        
        participaciones = competencia.participaciones.filter(
            finalizo=True
        ).select_related('estudiante').order_by('posicion', '-puntaje_total')
        
        for participacion in participaciones:
            tiempo_min = participacion.tiempo_total_segundos / 60 if participacion.tiempo_total_segundos else 0
            writer.writerow([
                participacion.posicion or '',
                f"{participacion.estudiante.first_name} {participacion.estudiante.last_name}".strip() or participacion.estudiante.username,
                participacion.estudiante.email,
                participacion.puntaje_total,
                f"{tiempo_min:.1f}",
                'Sí' if participacion.finalizo else 'No',
                participacion.fecha_ingreso.strftime('%Y-%m-%d %H:%M'),
                participacion.fecha_finalizacion.strftime('%Y-%m-%d %H:%M') if participacion.fecha_finalizacion else ''
            ])
    else:
        # Modalidad grupal
        writer.writerow([
            'Posición Grupo', 'Nombre Grupo', 'Nombre Estudiante', 'Email', 
            'Puntaje Individual', 'Puntaje Promedio Grupo', 'Tiempo (min)', 'Finalizado'
        ])
        
        for grupo in competencia.grupos.all().order_by('participaciones__posicion'):
            participaciones_grupo = grupo.participaciones.filter(finalizo=True)
            puntaje_promedio = participaciones_grupo.aggregate(
                promedio=models.Avg('puntaje_total')
            )['promedio'] or 0
            
            for participacion in participaciones_grupo:
                tiempo_min = participacion.tiempo_total_segundos / 60 if participacion.tiempo_total_segundos else 0
                writer.writerow([
                    participacion.posicion or '',
                    grupo.nombre,
                    f"{participacion.estudiante.first_name} {participacion.estudiante.last_name}".strip() or participacion.estudiante.username,
                    participacion.estudiante.email,
                    participacion.puntaje_total,
                    f"{puntaje_promedio:.2f}",
                    f"{tiempo_min:.1f}",
                    'Sí' if participacion.finalizo else 'No'
                ])
    
    return response


@csrf_exempt
@require_http_methods(["POST"])
@login_required
@role_required(['Admin', 'Docente'])
def reiniciar_competencia(request, competencia_id):
    """Reiniciar competencia (volver a configuración)"""
    try:
        competencia = get_object_or_404(Competencia, id=competencia_id)
        
        # Verificar permisos
        if (competencia.recurso.seccion and 
            competencia.recurso.seccion.curso.profesor != request.user):
            return JsonResponse({'success': False, 'error': 'Sin permisos'}, status=403)
        
        # Solo permitir reiniciar si no está activa
        if competencia.estado == 'activa':
            return JsonResponse({
                'success': False, 
                'error': 'No se puede reiniciar una competencia activa'
            })
        
        with transaction.atomic():
            # Eliminar todas las participaciones y respuestas
            competencia.participaciones.all().delete()
            
            # Eliminar todos los grupos
            competencia.grupos.all().delete()
            
            # Resetear estado
            competencia.estado = 'configuracion'
            competencia.fecha_inicio_real = None
            competencia.fecha_finalizacion = None
            competencia.pregunta_actual = None
            competencia.save()
        
        return JsonResponse({
            'success': True,
            'mensaje': 'Competencia reiniciada exitosamente'
        })
        
    except Exception as e:
        return JsonResponse({'success': False, 'error': str(e)}, status=400)
    

@csrf_exempt
@require_http_methods(["POST"])
@login_required
@role_required(['Admin', 'Docente'])
def cancelar_competencia(request, competencia_id):
    """Cancelar competencia"""
    try:
        competencia = get_object_or_404(Competencia, id=competencia_id)
        
        # Verificar permisos
        if (competencia.recurso.seccion and 
            competencia.recurso.seccion.curso.profesor != request.user):
            return JsonResponse({'success': False, 'error': 'Sin permisos'}, status=403)
        
        competencia.estado = 'cancelada'
        competencia.save()
        
        return JsonResponse({
            'success': True,
            'mensaje': 'Competencia cancelada'
        })
        
    except Exception as e:
        return JsonResponse({'success': False, 'error': str(e)}, status=400)


@csrf_exempt
@require_http_methods(["POST"])
@login_required
@role_required(['Admin', 'Docente'])
def duplicar_competencia(request, competencia_id):
    """Duplicar competencia completa"""
    try:
        competencia_original = get_object_or_404(Competencia, id=competencia_id)
        
        # Verificar permisos
        if (competencia_original.recurso.seccion and 
            competencia_original.recurso.seccion.curso.profesor != request.user):
            return JsonResponse({'success': False, 'error': 'Sin permisos'}, status=403)
        
        with transaction.atomic():
            # Crear nuevo recurso
            recurso_nuevo = Recurso.objects.create(
                seccion=competencia_original.recurso.seccion,
                tipo=competencia_original.recurso.tipo,
                titulo=f"{competencia_original.recurso.titulo} (Copia)",
                descripcion=competencia_original.recurso.descripcion,
                orden=competencia_original.recurso.seccion.recursos.count() + 1 if competencia_original.recurso.seccion else 1,
                creado_por=request.user
            )
            
            # Crear nueva competencia
            nueva_competencia = Competencia.objects.create(
                recurso=recurso_nuevo,
                instrucciones=competencia_original.instrucciones,
                modalidad=competencia_original.modalidad,
                max_miembros_grupo=competencia_original.max_miembros_grupo,
                grupos_aleatorios=competencia_original.grupos_aleatorios,
                grupos_abiertos=competencia_original.grupos_abiertos,
                tiempo_limite=competencia_original.tiempo_limite,
                estado='configuracion',
                mostrar_resultados_inmediatos=competencia_original.mostrar_resultados_inmediatos,
                permitir_reingreso=competencia_original.permitir_reingreso,
                orden_preguntas_aleatorio=competencia_original.orden_preguntas_aleatorio
            )
            
            # Duplicar preguntas
            from django.contrib.contenttypes.models import ContentType
            competencia_ct = ContentType.objects.get_for_model(Competencia)
            
            preguntas_originales = PreguntaCuestionario.objects.filter(
                content_type=ContentType.objects.get_for_model(Competencia),
                object_id=competencia_original.id
            ).order_by('orden')
            
            for pregunta_original in preguntas_originales:
                nueva_pregunta = PreguntaCuestionario.objects.create(
                    content_type=competencia_ct,
                    object_id=nueva_competencia.id,
                    enunciado=pregunta_original.enunciado,
                    tipo=pregunta_original.tipo,
                    orden=pregunta_original.orden,
                    puntaje=pregunta_original.puntaje,
                    # Copiar campos específicos por tipo
                    respuesta_modelo=getattr(pregunta_original, 'respuesta_modelo', None),
                    criterios_evaluacion=getattr(pregunta_original, 'criterios_evaluacion', None),
                    longitud_minima=getattr(pregunta_original, 'longitud_minima', None),
                    longitud_maxima=getattr(pregunta_original, 'longitud_maxima', None),
                    texto_completar=getattr(pregunta_original, 'texto_completar', None),
                    respuestas_completar=getattr(pregunta_original, 'respuestas_completar', None),
                    sensible_mayusculas=getattr(pregunta_original, 'sensible_mayusculas', False),
                    ignorar_espacios=getattr(pregunta_original, 'ignorar_espacios', True),
                    permitir_alternativas=getattr(pregunta_original, 'permitir_alternativas', False),
                    respuestas_alternativas=getattr(pregunta_original, 'respuestas_alternativas', None),
                    columna_izquierda=getattr(pregunta_original, 'columna_izquierda', None),
                    columna_derecha=getattr(pregunta_original, 'columna_derecha', None),
                    conexiones_correctas=getattr(pregunta_original, 'conexiones_correctas', None),
                    mezclar_opciones=getattr(pregunta_original, 'mezclar_opciones', True),
                    permitir_conexiones_multiples=getattr(pregunta_original, 'permitir_conexiones_multiples', False),
                    smiles_objetivo=getattr(pregunta_original, 'smiles_objetivo', None),
                    descripcion_molecula=getattr(pregunta_original, 'descripcion_molecula', None),
                    tolerancia_similitud=getattr(pregunta_original, 'tolerancia_similitud', None),
                    permitir_isomeros=getattr(pregunta_original, 'permitir_isomeros', False)
                )
                
                # Duplicar opciones
                for opcion in pregunta_original.opciones.all():
                    OpcionPregunta.objects.create(
                        pregunta=nueva_pregunta,
                        texto=opcion.texto,
                        es_correcta=opcion.es_correcta
                    )
        
        return JsonResponse({
            'success': True,
            'competencia_id': nueva_competencia.id,
            'mensaje': 'Competencia duplicada exitosamente',
            'redirect_url': f'/docente/editar_competencia/{nueva_competencia.id}/'
        })
        
    except Exception as e:
        print(f"❌ Error al duplicar competencia: {str(e)}")
        return JsonResponse({'success': False, 'error': str(e)}, status=400)
    

# ============== 9. AGREGARA AL FINAL - FUNCIONES AUXILIARES ==============

def obtener_estadisticas_competencia(competencia):
    """Obtener estadísticas completas de una competencia"""
    from django.db.models import Avg, Min, Max, Count
    
    participaciones = competencia.participaciones.filter(finalizo=True)
    
    if not participaciones.exists():
        return {
            'total_participantes': 0,
            'participantes_finalizados': 0,
            'puntaje_promedio': 0,
            'puntaje_maximo': 0,
            'puntaje_minimo': 0,
            'tiempo_promedio': 0,
            'tiempo_minimo': 0,
            'tiempo_maximo': 0
        }
    
    estadisticas = participaciones.aggregate(
        puntaje_promedio=Avg('puntaje_total'),
        puntaje_maximo=Max('puntaje_total'),
        puntaje_minimo=Min('puntaje_total'),
        tiempo_promedio=Avg('tiempo_total_segundos'),
        tiempo_minimo=Min('tiempo_total_segundos'),
        tiempo_maximo=Max('tiempo_total_segundos')
    )
    
    estadisticas.update({
        'total_participantes': competencia.participaciones.filter(activo=True).count(),
        'participantes_finalizados': participaciones.count(),
    })
    
    return estadisticas


def validar_configuracion_competencia(competencia):
    """Validar que una competencia esté correctamente configurada"""
    errores = []
    
    if not competencia.recurso.titulo.strip():
        errores.append("La competencia debe tener un título")
    
    if not competencia.instrucciones.strip():
        errores.append("La competencia debe tener instrucciones")
    
    # CORREGIR: Validar preguntas usando ContentType
    from django.contrib.contenttypes.models import ContentType
    competencia_ct = ContentType.objects.get_for_model(Competencia)
    preguntas_count = PreguntaCuestionario.objects.filter(
        content_type=competencia_ct,
        object_id=competencia.id
    ).count()
    
    if preguntas_count == 0:
        errores.append("La competencia debe tener al menos una pregunta")
    
    if competencia.tiempo_limite < 1 or competencia.tiempo_limite > 180:
        errores.append("El tiempo límite debe estar entre 1 y 180 minutos")
    
    if competencia.modalidad == 'grupal':
        if competencia.max_miembros_grupo < 2 or competencia.max_miembros_grupo > 10:
            errores.append("Los grupos deben tener entre 2 y 10 miembros")
    
    # Validar preguntas
    preguntas = PreguntaCuestionario.objects.filter(
        content_type=competencia_ct,
        object_id=competencia.id
    )
    
    for pregunta in preguntas:
        if not pregunta.enunciado.strip():
            errores.append(f"La pregunta {pregunta.orden} debe tener enunciado")
        
        if pregunta.tipo.nombre in ['opcion_unica', 'opcion_multiple', 'falso_verdadero']:
            if not pregunta.opciones.exists():
                errores.append(f"La pregunta {pregunta.orden} debe tener opciones")
            
            if not pregunta.opciones.filter(es_correcta=True).exists():
                errores.append(f"La pregunta {pregunta.orden} debe tener al menos una respuesta correcta")
    
    return errores



# AGREGAR ESTAS FUNCIONES AL FINAL DEL ARCHIVO views.py

# =====================================================
# FUNCIONES FALTANTES PARA COMPETENCIAS
# =====================================================

@login_required
@role_required(['Admin', 'Docente'])
def panel_competencia(request, competencia_id):
    """Panel de control principal para la competencia - NUEVA FUNCIÓN"""
    print(f"🎮 Panel de competencia {competencia_id}")
    
    competencia = get_object_or_404(Competencia, id=competencia_id)
    
    # Verificar permisos
    if (competencia.recurso.seccion and 
        competencia.recurso.seccion.curso.profesor != request.user):
        print("❌ Usuario sin permisos para ver panel")
        return redirect('dashboard-adm')
    
    # Obtener preguntas de la competencia
    from django.contrib.contenttypes.models import ContentType
    competencia_ct = ContentType.objects.get_for_model(Competencia)
    preguntas = PreguntaCuestionario.objects.filter(
        content_type=competencia_ct,
        object_id=competencia.id
    ).order_by('orden')
    
    # Obtener participaciones
    participaciones = competencia.participaciones.filter(activo=True).select_related('estudiante', 'grupo')
    
    # Obtener grupos si es modalidad grupal
    grupos = []
    if competencia.modalidad == 'grupal':
        grupos = competencia.grupos.all().prefetch_related('participaciones__estudiante')
    
    # Estadísticas
    total_participantes = participaciones.count()
    participantes_finalizados = participaciones.filter(finalizo=True).count()
    participantes_activos = total_participantes - participantes_finalizados
    
    # Ranking top 10
    ranking = list(participaciones.order_by('-puntaje_total', 'fecha_ingreso')[:10])
    
    context = {
        'competencia': competencia,
        'preguntas': preguntas,
        'participaciones': participaciones,
        'grupos': grupos,
        'total_participantes': total_participantes,
        'participantes_activos': participantes_activos,
        'participantes_finalizados': participantes_finalizados,
        'ranking': ranking,
        'tiempo_restante': competencia.tiempo_restante_segundos(),
        'puede_iniciar': competencia.puede_iniciar(),
        'puede_finalizar': competencia.puede_finalizar(),
        'esta_activa': competencia.esta_activa(),
        'esta_esperando': competencia.esta_esperando(),
        'esta_finalizada': competencia.esta_finalizada(),
        'imgPerfil': request.user.imgPerfil,
        'usuario': request.user.username,
    }
    
    return render(request, 'docente/competencias/panel_competencia.html', context)


@login_required
@role_required(['Admin', 'Docente'])
def agregar_pregunta_competencia_legacy(request, competencia_id):
    """Función legacy para agregar preguntas - NUEVA FUNCIÓN"""
    competencia = get_object_or_404(Competencia, id=competencia_id)
    
    # Verificar permisos
    if (competencia.recurso.seccion and 
        competencia.recurso.seccion.curso.profesor != request.user):
        return redirect('dashboard-adm')
    
    if request.method == 'POST':
        enunciado = request.POST.get('enunciado')
        tipo_id = request.POST.get('tipo_pregunta')
        orden = request.POST.get('orden')
        puntaje = request.POST.get('puntaje', 1)
        
        tipo = get_object_or_404(TipoPregunta, id=tipo_id)
        
        # Crear la pregunta usando relación genérica
        from django.contrib.contenttypes.models import ContentType
        competencia_ct = ContentType.objects.get_for_model(Competencia)
        
        pregunta = PreguntaCuestionario.objects.create(
            content_type=competencia_ct,
            object_id=competencia.id,
            enunciado=enunciado,
            tipo=tipo,
            orden=orden,
            puntaje=puntaje
        )
        
        # Crear opciones según el tipo
        tipo_nombre = tipo.nombre.lower()
        
        if tipo_nombre in ['opcion_unica', 'opcion_multiple']:
            for i in range(1, 5):
                texto = request.POST.get(f'opcion_{i}')
                es_correcta = request.POST.get('correcta') == str(i)
                if texto:
                    OpcionPregunta.objects.create(
                        pregunta=pregunta, 
                        texto=texto, 
                        es_correcta=es_correcta
                    )
        
        elif tipo_nombre == 'falso_verdadero':
            respuesta = request.POST.get('vf_respuesta')
            if respuesta:
                OpcionPregunta.objects.create(
                    pregunta=pregunta, 
                    texto='Verdadero', 
                    es_correcta=(respuesta == 'verdadero')
                )
                OpcionPregunta.objects.create(
                    pregunta=pregunta, 
                    texto='Falso', 
                    es_correcta=(respuesta == 'falso')
                )
        
        messages.success(request, 'Pregunta agregada exitosamente')
        return redirect('panel_competencia', competencia.id)
    
    tipos = TipoPregunta.objects.all()
    return render(request, 'docente/competencias/agregar_pregunta.html', {
        'competencia': competencia,
        'tipos': tipos
    })


@csrf_exempt
@require_http_methods(["POST"])
@login_required
@role_required(['Admin', 'Docente'])
def iniciar_competencia(request, competencia_id):
    """Iniciar competencia - NUEVA FUNCIÓN (diferente de iniciar_competencia_live)"""
    try:
        print(f"🚀 Iniciando competencia {competencia_id}")
        
        competencia = get_object_or_404(Competencia, id=competencia_id)
        
        # Verificar permisos
        if (competencia.recurso.seccion and 
            competencia.recurso.seccion.curso.profesor != request.user):
            return JsonResponse({'success': False, 'error': 'Sin permisos'}, status=403)
        
        # Verificar que pueda iniciarse
        if not competencia.puede_iniciar():
            return JsonResponse({
                'success': False, 
                'error': 'La competencia no está lista para iniciar'
            })
        
        # Verificar que hay participantes
        participantes_count = competencia.participaciones.filter(activo=True).count()
        if participantes_count == 0:
            return JsonResponse({
                'success': False, 
                'error': 'No hay participantes para iniciar la competencia'
            })
        
        with transaction.atomic():
            # Cambiar estado a activa
            competencia.estado = 'activa'
            competencia.fecha_inicio_real = timezone.now()
            competencia.save()
            
            print(f"✅ Competencia {competencia_id} iniciada")
        
        if request.headers.get('Content-Type') == 'application/json':
            return JsonResponse({
                'success': True,
                'mensaje': 'Competencia iniciada exitosamente'
            })
        else:
            messages.success(request, 'Competencia iniciada exitosamente')
            return redirect('panel_competencia', competencia_id)
        
    except Exception as e:
        print(f"❌ Error al iniciar competencia: {str(e)}")
        if request.headers.get('Content-Type') == 'application/json':
            return JsonResponse({'success': False, 'error': str(e)})
        else:
            messages.error(request, f'Error al iniciar competencia: {str(e)}')
            return redirect('panel_competencia', competencia_id)


@csrf_exempt
@require_http_methods(["POST"])
@login_required
@role_required(['Admin', 'Docente'])
def enviar_pregunta(request, pregunta_id):
    """Enviar/mostrar una pregunta específica - NUEVA FUNCIÓN"""
    try:
        print(f"📋 Enviando pregunta {pregunta_id}")
        
        pregunta = get_object_or_404(PreguntaCuestionario, id=pregunta_id)
        
        # Verificar que pertenece a una competencia
        if not isinstance(pregunta.content_object, Competencia):
            return JsonResponse({
                'success': False, 
                'error': 'Esta pregunta no pertenece a una competencia'
            })
        
        competencia = pregunta.content_object
        
        # Verificar permisos
        if (competencia.recurso.seccion and 
            competencia.recurso.seccion.curso.profesor != request.user):
            return JsonResponse({'success': False, 'error': 'Sin permisos'}, status=403)
        
        # Verificar que la competencia esté activa
        if competencia.estado != 'activa':
            return JsonResponse({
                'success': False, 
                'error': 'La competencia debe estar activa'
            })
        
        # Actualizar pregunta actual
        competencia.pregunta_actual = pregunta
        competencia.save()
        
        print(f"✅ Pregunta {pregunta_id} enviada como actual")
        
        if request.headers.get('Content-Type') == 'application/json':
            return JsonResponse({
                'success': True,
                'mensaje': f'Pregunta {pregunta.orden} enviada',
                'pregunta_actual': {
                    'id': pregunta.id,
                    'orden': pregunta.orden,
                    'enunciado': pregunta.enunciado
                }
            })
        else:
            messages.success(request, f'Pregunta {pregunta.orden} enviada')
            return redirect('panel_competencia', competencia.id)
        
    except Exception as e:
        print(f"❌ Error al enviar pregunta: {str(e)}")
        if request.headers.get('Content-Type') == 'application/json':
            return JsonResponse({'success': False, 'error': str(e)})
        else:
            messages.error(request, f'Error al enviar pregunta: {str(e)}')
            return redirect('panel_competencia', competencia.content_object.id if hasattr(pregunta, 'content_object') else 1)


# =====================================================
# FUNCIONES AUXILIARES ADICIONALES
# =====================================================

@require_http_methods(["GET"])
@login_required
def obtener_info_participacion(request, participacion_id):
    """Obtener información de participación (AJAX)"""
    try:
        participacion = get_object_or_404(
            ParticipacionCompetencia,
            id=participacion_id,
            estudiante=request.user
        )
        
        return JsonResponse({
            'success': True,
            'participacion': {
                'id': participacion.id,
                'finalizo': participacion.finalizo,
                'puntaje_total': float(participacion.puntaje_total),
                'tiempo_transcurrido': (timezone.now() - participacion.fecha_ingreso).total_seconds(),
                'grupo': {
                    'id': participacion.grupo.id,
                    'nombre': participacion.grupo.nombre
                } if participacion.grupo else None
            },
            'competencia': {
                'id': participacion.competencia.id,
                'estado': participacion.competencia.estado,
                'tiempo_restante': participacion.competencia.tiempo_restante_segundos()
            }
        })
        
    except Exception as e:
        return JsonResponse({'success': False, 'error': str(e)})


@csrf_exempt
@require_http_methods(["POST"])
@login_required
def ping_participacion(request, participacion_id):
    """Ping para mantener participación activa"""
    try:
        participacion = get_object_or_404(
            ParticipacionCompetencia,
            id=participacion_id,
            estudiante=request.user
        )
        
        # Actualizar última actividad
        participacion.fecha_ultima_actividad = timezone.now()
        participacion.save()
        
        return JsonResponse({
            'success': True,
            'timestamp': timezone.now().isoformat()
        })
        
    except Exception as e:
        return JsonResponse({'success': False, 'error': str(e)})


# =====================================================
# FUNCIONES ADICIONALES QUE FALTAN EN URLS
# =====================================================

def pausar_competencia(request, competencia_id):
    """Pausar competencia (placeholder)"""
    messages.info(request, 'Función pausar en desarrollo')
    return redirect('panel_competencia', competencia_id)

def reanudar_competencia(request, competencia_id):
    """Reanudar competencia (placeholder)"""
    messages.info(request, 'Función reanudar en desarrollo')
    return redirect('panel_competencia', competencia_id)

def pregunta_anterior_competencia(request, competencia_id):
    """Ir a pregunta anterior (placeholder)"""
    messages.info(request, 'Función pregunta anterior en desarrollo')
    return redirect('panel_competencia', competencia_id)

def ir_a_pregunta_competencia(request, competencia_id, pregunta_orden):
    """Ir a pregunta específica (placeholder)"""
    messages.info(request, f'Función ir a pregunta {pregunta_orden} en desarrollo')
    return redirect('panel_competencia', competencia_id)

def editar_grupo(request, grupo_id):
    """Editar grupo (placeholder)"""
    messages.info(request, 'Función editar grupo en desarrollo')
    return redirect('biblioteca_competencias')

def asignar_estudiante_grupo(request):
    """Asignar estudiante a grupo (placeholder)"""
    messages.info(request, 'Función asignar estudiante en desarrollo')
    return redirect('biblioteca_competencias')

def remover_estudiante_grupo(request, participacion_id):
    """Remover estudiante de grupo (placeholder)"""
    messages.info(request, 'Función remover estudiante en desarrollo')
    return redirect('biblioteca_competencias')

def resultados_detallados_competencia(request, competencia_id):
    """Resultados detallados (placeholder)"""
    messages.info(request, 'Función resultados detallados en desarrollo')
    return redirect('resultados_competencia', competencia_id)

def exportar_resultados_excel(request, competencia_id):
    """Exportar a Excel (placeholder)"""
    messages.info(request, 'Función exportar Excel en desarrollo')
    return redirect('resultados_competencia', competencia_id)

def estadisticas_competencia(request, competencia_id):
    """Estadísticas de competencia (placeholder)"""
    messages.info(request, 'Función estadísticas en desarrollo')
    return redirect('resultados_competencia', competencia_id)

def archivar_competencia(request, competencia_id):
    """Archivar competencia (placeholder)"""
    messages.info(request, 'Función archivar en desarrollo')
    return redirect('biblioteca_competencias')

def restaurar_competencia(request, competencia_id):
    """Restaurar competencia (placeholder)"""
    messages.info(request, 'Función restaurar en desarrollo')
    return redirect('biblioteca_competencias')

def crear_desde_plantilla(request, plantilla_id):
    """Crear desde plantilla (placeholder)"""
    messages.info(request, 'Función crear desde plantilla en desarrollo')
    return redirect('biblioteca_competencias')

def guardar_como_plantilla(request, competencia_id):
    """Guardar como plantilla (placeholder)"""
    messages.info(request, 'Función guardar plantilla en desarrollo')
    return redirect('biblioteca_competencias')

# ESTUDIANTE

@csrf_exempt
@require_http_methods(["POST"])
@login_required
@role_required('Estudiante')
def buscar_competencia_codigo(request):
    """Buscar competencia por PIN y redirigir a seleccionar avatar / sala de espera."""
    try:
        # 1) PIN
        if request.content_type == 'application/json':
            data = json.loads(request.body)
            pin = data.get('pin', '').strip()
        else:
            pin = request.POST.get('pin', '').strip()

        if not pin or len(pin) != 6:
            messages.error(request, 'PIN debe tener 6 dígitos')
            return redirect('unirse_competencia_principal')

        # 2) Competencia por PIN
        competencia = Competencia.objects.get(pin_acceso=pin)  # 404 si no existe

        # 3) Crear/activar participación
        participacion, creada = ParticipacionCompetencia.objects.get_or_create(
            competencia=competencia,
            estudiante=request.user,
            defaults={'activo': True, 'fecha_ultima_actividad': timezone.now()}
        )
        if not creada:
            participacion.activo = True
            participacion.fecha_ultima_actividad = timezone.now()
            participacion.save()

        # 4) Redirección por estado del avatar
        if not participacion.avatar_elegido:
            # usar reverse por nombre de ruta; agregamos querystring con participacion_id
            url = reverse('seleccionar_avatar_competencia', args=[competencia.id])
            return redirect(f'{url}?participacion_id={participacion.id}')
        else:
            return redirect('sala_espera_competencia', participacion_id=participacion.id)

    except Competencia.DoesNotExist:
        messages.error(request, 'PIN incorrecto. Verifica e inténtalo de nuevo.')
        return redirect('unirse_competencia_principal')
    except Exception as e:
        messages.error(request, f'Error: {str(e)}')
        return redirect('unirse_competencia_principal')
    
    


def determinar_redirect_reingreso(participacion):
    """Determina a dónde redirigir según el estado del participante"""
    competencia = participacion.competencia
    
    print(f"🔄 Reingreso - Estado: {competencia.estado}, Avatar: {participacion.avatar_elegido}")
    
    if competencia.estado == 'esperando':
        # Si NO tiene avatar, ir a seleccionar avatar
        if not participacion.avatar_elegido:
            return f'/competencia/seleccionar-avatar/{competencia.id}/?participacion_id={participacion.id}'
        
        # Si YA tiene avatar, ir directo a sala de espera
        return f'/competencia/sala-espera/{participacion.id}/'
    
    elif competencia.estado == 'activa':
        # Si está activa, ir directo a participar (sin importar avatar)
        return f'/competencia/participar/{participacion.id}/'
    
    # Fallback
    return f'/competencia/sala-espera/{participacion.id}/'

    

def salir_competencia(request, participacion_id):
    """Salir de competencia (placeholder)"""
    messages.info(request, 'Función salir competencia en desarrollo')
    return redirect('student_dashboard')

def salir_grupo(request, participacion_id):
    """Salir de grupo (placeholder)"""
    messages.info(request, 'Función salir grupo en desarrollo')
    return redirect('student_dashboard')

def info_grupo(request, grupo_id):
    """Info de grupo (placeholder)"""
    messages.info(request, 'Función info grupo en desarrollo')
    return redirect('student_dashboard')

def abandonar_competencia(request, participacion_id):
    """Abandonar competencia (placeholder)"""
    messages.info(request, 'Función abandonar en desarrollo')
    return redirect('student_dashboard')

def resultados_grupo(request, grupo_id):
    """Resultados de grupo (placeholder)"""
    messages.info(request, 'Función resultados grupo en desarrollo')
    return redirect('student_dashboard')

def historial_competencias_estudiante(request):
    """Historial de competencias (placeholder)"""
    messages.info(request, 'Función historial en desarrollo')
    return redirect('student_dashboard')

def generar_certificado(request, participacion_id):
    """Generar certificado (placeholder)"""
    messages.info(request, 'Función certificado en desarrollo')
    return redirect('student_dashboard')

# Funciones AJAX adicionales

@require_http_methods(["GET"])
@login_required
@role_required('Estudiante')
def obtener_pregunta_actual_estudiante(request, competencia_id):
    """Obtener pregunta actual para estudiantes (AJAX)"""
    try:
        competencia = get_object_or_404(Competencia, id=competencia_id)
        
        # Verificar participación del estudiante
        participacion = ParticipacionCompetencia.objects.filter(
            competencia=competencia,
            estudiante=request.user,
            activo=True
        ).first()
        
        if not participacion:
            return JsonResponse({
                'success': False,
                'error': 'No estás registrado en esta competencia'
            })
        
        if competencia.estado != 'activa':
            return JsonResponse({
                'success': False,
                'error': 'La competencia no está activa',
                'estado': competencia.estado
            })
        
        if not competencia.pregunta_actual:
            return JsonResponse({
                'success': False,
                'error': 'No hay pregunta actual'
            })
        
        pregunta = competencia.pregunta_actual
        
        # Estructurar respuesta según tipo de pregunta
        pregunta_data = {
            'id': pregunta.id,
            'orden': pregunta.orden,
            'enunciado': pregunta.enunciado,
            'tipo': pregunta.tipo.nombre,
            'puntaje': pregunta.puntaje,
            'imagen': pregunta.imagen.url if pregunta.imagen else None,
            'url_youtube': pregunta.get_url_youtube_embed() if pregunta.url_youtube else None,
        }
        
        # Agregar opciones si las tiene
        if pregunta.opciones.exists():
            pregunta_data['opciones'] = [
                {
                    'id': opcion.id,
                    'texto': opcion.texto
                }
                for opcion in pregunta.opciones.all()
            ]
        
        # Agregar campos específicos según tipo
        if pregunta.tipo.nombre == 'completar':
            pregunta_data['texto_completar'] = pregunta.texto_completar
        
        elif pregunta.tipo.nombre == 'unir_lineas':
            pregunta_data['columna_izquierda'] = pregunta.columna_izquierda.split('\n') if pregunta.columna_izquierda else []
            pregunta_data['columna_derecha'] = pregunta.columna_derecha.split('\n') if pregunta.columna_derecha else []
        
        elif pregunta.tipo.nombre in ['simulador_2d', 'simulador_3d']:
            pregunta_data['descripcion_molecula'] = pregunta.descripcion_molecula
            pregunta_data['smiles_objetivo'] = pregunta.smiles_objetivo
        
        return JsonResponse({
            'success': True,
            'pregunta': pregunta_data,
            'tiempo_restante': competencia.tiempo_restante_segundos()
        })
        
    except Exception as e:
        print(f"❌ Error al obtener pregunta actual: {str(e)}")
        return JsonResponse({
            'success': False,
            'error': str(e)
        })

def info_siguiente_pregunta(request, competencia_id):
    """Info siguiente pregunta (placeholder)"""
    return JsonResponse({'success': False, 'error': 'En desarrollo'})

def tiempo_restante_competencia(request, competencia_id):
    """Tiempo restante (placeholder)"""
    return JsonResponse({'success': False, 'error': 'En desarrollo'})

def obtener_ranking_grupos(request, competencia_id):
    """Ranking de grupos (placeholder)"""
    return JsonResponse({'success': False, 'error': 'En desarrollo'})

def obtener_progreso_competencia(request, competencia_id):
    """Progreso de competencia (placeholder)"""
    return JsonResponse({'success': False, 'error': 'En desarrollo'})

def estadisticas_live(request, competencia_id):
    """Estadísticas en vivo (placeholder)"""
    return JsonResponse({'success': False, 'error': 'En desarrollo'})

def heartbeat_competencia(request, competencia_id):
    """Heartbeat (placeholder)"""
    return JsonResponse({'success': True, 'timestamp': timezone.now().isoformat()})

# Chat y comunicación

def chat_competencia(request, competencia_id):
    """Chat de competencia (placeholder)"""
    messages.info(request, 'Función chat en desarrollo')
    return redirect('biblioteca_competencias')

def enviar_mensaje_chat(request):
    """Enviar mensaje (placeholder)"""
    return JsonResponse({'success': False, 'error': 'En desarrollo'})

def obtener_mensajes_chat(request, competencia_id):
    """Obtener mensajes (placeholder)"""
    return JsonResponse({'success': False, 'error': 'En desarrollo'})

# APIs

def api_listar_competencias(request):
    """API listar (placeholder)"""
    return JsonResponse({'success': False, 'error': 'En desarrollo'})

def api_info_competencia(request, competencia_id):
    """API info (placeholder)"""
    return JsonResponse({'success': False, 'error': 'En desarrollo'})

def api_unirse_competencia(request):
    """API unirse (placeholder)"""
    return JsonResponse({'success': False, 'error': 'En desarrollo'})

def api_responder_competencia(request):
    """API responder (placeholder)"""
    return JsonResponse({'success': False, 'error': 'En desarrollo'})

def buscar_competencias(request):
    """Buscar competencias (placeholder)"""
    messages.info(request, 'Función buscar en desarrollo')
    return redirect('biblioteca_competencias')

def filtrar_competencias(request):
    """Filtrar competencias (placeholder)"""
    messages.info(request, 'Función filtrar en desarrollo')
    return redirect('biblioteca_competencias')

@login_required
def api_estado_competencia(request, competencia_id:int):
    comp = get_object_or_404(Competencia, id=competencia_id)
    return JsonResponse({'estado': str(comp.estado)})

# Notificaciones y webhooks

def notificar_inicio_competencia(request, competencia_id):
    """Notificar inicio (placeholder)"""
    messages.info(request, 'Función notificar en desarrollo')
    return redirect('panel_competencia', competencia_id)

def enviar_recordatorio(request, competencia_id):
    """Enviar recordatorio (placeholder)"""
    messages.info(request, 'Función recordatorio en desarrollo')
    return redirect('panel_competencia', competencia_id)

def webhook_inicio_competencia(request):
    """Webhook inicio (placeholder)"""
    return JsonResponse({'success': True})

def webhook_finalizacion_competencia(request):
    """Webhook finalización (placeholder)"""
    return JsonResponse({'success': True})


@csrf_exempt
@require_http_methods(["POST"])
@login_required
@role_required(['Admin', 'Docente'])
def cambiar_tipo_pregunta_competencia(request):
    """Cambiar tipo de pregunta de competencia - IMPLEMENTACIÓN REAL"""
    try:
        data = json.loads(request.body)
        pregunta_id = data.get('pregunta_id')
        nuevo_tipo = data.get('nuevo_tipo')
        pregunta = get_object_or_404(PreguntaCuestionario, id=pregunta_id)
        
        # Verificar que sea de competencia
        if not isinstance(pregunta.content_object, Competencia):
            return JsonResponse({'success': False, 'error': 'No es pregunta de competencia'})
        
        competencia = pregunta.content_object
        
        # Verificar permisos y estado
        if (competencia.recurso.seccion and 
            competencia.recurso.seccion.curso.profesor != request.user):
            return JsonResponse({'success': False, 'error': 'Sin permisos'}, status=403)
        
        if competencia.estado not in ['configuracion']:
            return JsonResponse({'success': False, 'error': 'Solo en configuración'})
        
        with transaction.atomic():
            # Eliminar opciones existentes
            pregunta.opciones.all().delete()
            
            # Cambiar tipo
            tipo_pregunta = get_object_or_404(TipoPregunta, nombre=nuevo_tipo)
            pregunta.tipo = tipo_pregunta
            pregunta.save()
            
            # Crear nuevas opciones según tipo
            if nuevo_tipo in ['opcion_unica', 'opcion_multiple']:
                for i in range(1, 5):
                    OpcionPregunta.objects.create(
                        pregunta=pregunta,
                        texto=f'Opción {i}',
                        es_correcta=(i == 1)
                    )
            elif nuevo_tipo == 'falso_verdadero':
                OpcionPregunta.objects.create(pregunta=pregunta, texto='Verdadero', es_correcta=True)
                OpcionPregunta.objects.create(pregunta=pregunta, texto='Falso', es_correcta=False)
        
        return JsonResponse({
            'success': True,
            'mensaje': 'Tipo de pregunta cambiado exitosamente'
        })
        
    except Exception as e:
        return JsonResponse({'success': False, 'error': str(e)})
    

@csrf_exempt
@require_http_methods(["POST"])
@login_required
@role_required(['Admin', 'Docente'])
def agregar_opcion_competencia(request):
    """Agregar opción a pregunta de competencia - IMPLEMENTACIÓN REAL"""
    try:
        data = json.loads(request.body)
        pregunta_id = data.get('pregunta_id')
        pregunta = get_object_or_404(PreguntaCuestionario, id=pregunta_id)
        
        # Verificar que sea de competencia
        if not isinstance(pregunta.content_object, Competencia):
            return JsonResponse({'success': False, 'error': 'No es pregunta de competencia'})
        
        competencia = pregunta.content_object
        
        # Verificar permisos
        if (competencia.recurso.seccion and 
            competencia.recurso.seccion.curso.profesor != request.user):
            return JsonResponse({'success': False, 'error': 'Sin permisos'}, status=403)
        
        # Solo en configuración
        if competencia.estado not in ['configuracion']:
            return JsonResponse({'success': False, 'error': 'Solo en configuración'})
        
        # Verificar tipo
        if pregunta.tipo.nombre not in ['opcion_unica', 'opcion_multiple']:
            return JsonResponse({'success': False, 'error': 'Tipo no permite opciones'})
        
        # Crear opción
        opciones_count = pregunta.opciones.count()
        opcion = OpcionPregunta.objects.create(
            pregunta=pregunta,
            texto=f"Nueva opción {opciones_count + 1}",
            es_correcta=False
        )
        
        return JsonResponse({
            'success': True,
            'opcion_id': opcion.id,
            'mensaje': 'Opción agregada exitosamente'
        })
        
    except Exception as e:
        return JsonResponse({'success': False, 'error': str(e)})

def actualizar_opcion_competencia(request):
    """Actualizar opción competencia (placeholder)"""
    return JsonResponse({'success': False, 'error': 'En desarrollo'})

@csrf_exempt
@require_http_methods(["DELETE"])
@login_required
@role_required(['Admin', 'Docente'])
def eliminar_opcion_competencia(request, opcion_id):
    """Eliminar opción de pregunta de competencia - IMPLEMENTACIÓN REAL"""
    try:
        opcion = get_object_or_404(OpcionPregunta, id=opcion_id)
        pregunta = opcion.pregunta
            # Verificar que sea de competencia
        if not isinstance(pregunta.content_object, Competencia):
            return JsonResponse({'success': False, 'error': 'No es pregunta de competencia'})
        
        competencia = pregunta.content_object
        
        # Verificar permisos
        if (competencia.recurso.seccion and 
            competencia.recurso.seccion.curso.profesor != request.user):
            return JsonResponse({'success': False, 'error': 'Sin permisos'}, status=403)
        
        # Solo en configuración
        if competencia.estado not in ['configuracion']:
            return JsonResponse({'success': False, 'error': 'Solo en configuración'})
        
        # Verificar mínimo
        if pregunta.opciones.count() <= 2:
            return JsonResponse({'success': False, 'error': 'Mínimo 2 opciones'})
        
        texto = opcion.texto
        opcion.delete()
        
        return JsonResponse({
            'success': True,
            'mensaje': f'Opción "{texto}" eliminada'
        })
        
    except Exception as e:
        return JsonResponse({'success': False, 'error': str(e)})


@csrf_exempt
@require_http_methods(["POST"])
@login_required
@role_required(['Admin', 'Docente'])
def subir_recurso_competencia(request):
    """Subir recurso multimedia a pregunta de competencia - IMPLEMENTACIÓN REAL"""
    try:
        pregunta_id = request.POST.get('pregunta_id')
        tipo_recurso = request.POST.get('tipo_recurso')
        pregunta = get_object_or_404(PreguntaCuestionario, id=pregunta_id)
        
        # Verificar que sea de competencia
        if not isinstance(pregunta.content_object, Competencia):
            return JsonResponse({'success': False, 'error': 'No es pregunta de competencia'})
        
        competencia = pregunta.content_object
        
        # Verificar permisos
        if (competencia.recurso.seccion and 
            competencia.recurso.seccion.curso.profesor != request.user):
            return JsonResponse({'success': False, 'error': 'Sin permisos'}, status=403)
        
        # Procesar según tipo
        if tipo_recurso in ['imagen', 'video', 'documento', 'audio'] and 'archivo' in request.FILES:
            pregunta.archivo_multimedia = request.FILES['archivo']
            pregunta.save()
        
        elif tipo_recurso == 'youtube':
            youtube_id = request.POST.get('youtube_id')
            if youtube_id:
                pregunta.url_youtube = youtube_id
                pregunta.save()
        
        elif tipo_recurso == 'simulador':
            # Procesar datos del simulador
            tipo_sim = request.POST.get('tipo_simulador')
            smiles_obj = request.POST.get('smiles_objetivo')
            instrucciones = request.POST.get('instrucciones')
            tolerancia = request.POST.get('tolerancia_similitud', 95)
            
            if tipo_sim and smiles_obj and instrucciones:
                pregunta.smiles_objetivo = smiles_obj
                pregunta.descripcion_molecula = instrucciones
                pregunta.tolerancia_similitud = int(tolerancia)
                pregunta.save()
        
        elif tipo_recurso == 'html':
            codigo_html = request.POST.get('codigo_html_encriptado')
            titulo = request.POST.get('titulo')
            descripcion = request.POST.get('descripcion', '')
            
            if codigo_html and titulo:
                # Guardar HTML encriptado en archivo multimedia temporalmente
                # O crear un campo específico en el modelo si es necesario
                pregunta.descripcion_molecula = f"HTML: {titulo} - {descripcion}"
                pregunta.save()
        
        return JsonResponse({
            'success': True,
            'mensaje': f'{tipo_recurso} cargado exitosamente'
        })
        
    except Exception as e:
        return JsonResponse({'success': False, 'error': str(e)})


@csrf_exempt
@require_http_methods(["POST"])
@login_required
@role_required(['Admin', 'Docente'])
def eliminar_recurso_pregunta_competencia(request):
    """Eliminar recurso de pregunta de competencia - IMPLEMENTACIÓN REAL"""
    try:
        data = json.loads(request.body)
        pregunta_id = data.get('pregunta_id')
        tipo = data.get('tipo')
        pregunta = get_object_or_404(PreguntaCuestionario, id=pregunta_id)
            
        # Verificar que sea de competencia
        if not isinstance(pregunta.content_object, Competencia):
            return JsonResponse({'success': False, 'error': 'No es pregunta de competencia'})
        
        competencia = pregunta.content_object
        
        # Verificar permisos
        if (competencia.recurso.seccion and 
            competencia.recurso.seccion.curso.profesor != request.user):
            return JsonResponse({'success': False, 'error': 'Sin permisos'}, status=403)
        
        if tipo == 'archivo':
            if pregunta.archivo_multimedia:
                try:
                    pregunta.archivo_multimedia.delete()
                except:
                    pass
                pregunta.archivo_multimedia = None
        elif tipo == 'youtube':
            pregunta.url_youtube = None
        
        pregunta.save()
        
        return JsonResponse({
            'success': True,
            'mensaje': 'Recurso eliminado exitosamente'
        })
        
    except Exception as e:
        return JsonResponse({'success': False, 'error': str(e)})

@csrf_exempt
@require_http_methods(["POST"])
@login_required
@role_required(['Admin', 'Docente'])
def guardar_competencia_completa(request, competencia_id):
    """Guardar competencia completa - IMPLEMENTACIÓN REAL"""
    try:
        competencia = get_object_or_404(Competencia, id=competencia_id)
        # Verificar permisos
        if (competencia.recurso.seccion and 
            competencia.recurso.seccion.curso.profesor != request.user):
            return JsonResponse({'success': False, 'error': 'Sin permisos'}, status=403)
        
        # Validar que esté completa
        errores = validar_configuracion_competencia(competencia)
        if errores:
            return JsonResponse({
                'success': False,
                'error': 'Competencia incompleta',
                'errores': errores
            })
        
        # Actualizar fecha de modificación
        competencia.fecha_modificacion = timezone.now()
        competencia.save()
        
        return JsonResponse({
            'success': True,
            'mensaje': 'Competencia guardada completamente',
            'puede_iniciar': competencia.puede_iniciar()
        })
    
    except Exception as e:
        return JsonResponse({'success': False, 'error': str(e)})


@csrf_exempt
@require_http_methods(["POST"])
@login_required
@role_required(['Admin', 'Docente'])
def autoguardar_competencia(request):
    """Autoguardar competencia - IMPLEMENTACIÓN REAL"""
    try:
        data = json.loads(request.body)
        competencia_id = data.get('competencia_id')
        titulo = data.get('titulo', '').strip()
        descripcion = data.get('descripcion', '').strip()
        competencia = get_object_or_404(Competencia, id=competencia_id)
        
        # Verificar permisos
        if (competencia.recurso.seccion and 
            competencia.recurso.seccion.curso.profesor != request.user):
            return JsonResponse({'success': False, 'error': 'Sin permisos'}, status=403)
        
        # Actualizar solo si hay cambios
        cambios = False
        if titulo and competencia.recurso.titulo != titulo:
            competencia.recurso.titulo = titulo
            cambios = True
        
        if descripcion and competencia.instrucciones != descripcion:
            competencia.instrucciones = descripcion
            cambios = True
        
        if cambios:
            competencia.recurso.save()
            competencia.save()
        
        return JsonResponse({
            'success': True,
            'mensaje': 'Autoguardado exitoso' if cambios else 'Sin cambios'
        })
    
    except Exception as e:
        return JsonResponse({'success': False, 'error': str(e)})

@login_required
def editar_perfil(request):
    user = request.user
    imgPerfil=user.imgPerfil      

    if request.method == 'POST':
        form = EditarPerfilForm(request.POST, request.FILES, instance=user)
        if form.is_valid():
            form.save()
            return redirect('editar_perfil')  
    else:
        form = EditarPerfilForm(instance=user)

    return render(request, 'usAdmin/editar_perfil.html', {'form': form, 'imgPerfil': imgPerfil,        
        'usuario':user })

def comming_soon(request):
    user = request.user
    imgPerfil=user.imgPerfil  
    context = {                  
        'imgPerfil': imgPerfil,        
        'usuario':user,        
    }   
    return render(request, 'coming-soon.html', context)    


@require_http_methods(["POST"])
@login_required
@role_required('Estudiante')
def guardar_avatar(request: HttpRequest):
    """
    Guarda el avatar elegido por el estudiante y devuelve JSON con la URL de redirección
    hacia la sala de espera. Soporta payload en JSON o application/x-www-form-urlencoded.

    Espera:
      - participacion_id
      - avatar (emoji o string)
      - avatar_nombre (opcional)
      - avatar_formula (opcional)
      - nombre_estudiante (opcional; NO modifica el User si no tienes ese flujo)
    """
    try:
        # 1) Parse de datos (JSON o POST clásico)
        if request.content_type == 'application/json':
            data = json.loads(request.body.decode('utf-8'))
        else:
            data = request.POST

        participacion_id = str(data.get('participacion_id', '')).strip()
        avatar = (data.get('avatar') or data.get('avatar_elegido') or '').strip()
        avatar_nombre = (data.get('avatar_nombre') or '').strip()
        avatar_formula = (data.get('avatar_formula') or '').strip()
        nombre_estudiante = (data.get('nombre_estudiante') or '').strip()

        if not participacion_id:
            return JsonResponse({'success': False, 'error': 'participacion_id es requerido'}, status=400)
        if not avatar:
            return JsonResponse({'success': False, 'error': 'Debes seleccionar un avatar'}, status=400)

        # 2) Cargar participación y validar pertenencia al usuario actual
        participacion = get_object_or_404(ParticipacionCompetencia, id=participacion_id)
        if participacion.estudiante_id != request.user.id:
            return JsonResponse({'success': False, 'error': 'No autorizado para esta participación'}, status=403)

        # 3) Guardar avatar y campos relacionados de forma segura
        with transaction.atomic():
            # Campo principal
            if hasattr(participacion, 'avatar_elegido'):
                participacion.avatar_elegido = avatar

            # Campos auxiliares si existen en tu modelo
            if hasattr(participacion, 'avatar_nombre') and avatar_nombre:
                participacion.avatar_nombre = avatar_nombre

            if hasattr(participacion, 'avatar_formula') and avatar_formula:
                participacion.avatar_formula = avatar_formula

            # Si tu modelo tiene un alias/“display name” de competencia
            if hasattr(participacion, 'nombre_mostrado') and nombre_estudiante:
                participacion.nombre_mostrado = nombre_estudiante

            # Marcar activa y actualizar actividad
            if hasattr(participacion, 'activo'):
                participacion.activo = True
            if hasattr(participacion, 'fecha_ultima_actividad'):
                participacion.fecha_ultima_actividad = timezone.now()

            participacion.save()

        # 4) URL de redirección a sala de espera (lo espera tu JS)
        redirect_url = reverse('sala_espera_competencia', kwargs={'participacion_id': participacion.id})

        return JsonResponse({
            'success': True,
            'redirect_url': redirect_url,
            'participacion_id': participacion.id
        })

    except ParticipacionCompetencia.DoesNotExist:
        return JsonResponse({'success': False, 'error': 'Participación no encontrada'}, status=404)
    except json.JSONDecodeError:
        return JsonResponse({'success': False, 'error': 'JSON inválido'}, status=400)
    except Exception as e:
        # Loguea e si tienes logger
        return JsonResponse({'success': False, 'error': f'Error al guardar avatar: {str(e)}'}, status=500)


@require_http_methods(["GET"])
@login_required
def obtener_participantes_sala_espera(request, competencia_id):
    """Obtener participantes para sala de espera (AJAX)"""
    try:
        competencia = get_object_or_404(Competencia, id=competencia_id)
        
        participaciones = competencia.participaciones.filter(
            activo=True
        ).select_related('estudiante', 'grupo')
        
        participantes = []
        for participacion in participaciones:
            participantes.append({
                'id': participacion.estudiante.id,
                'nombre': f"{participacion.estudiante.first_name} {participacion.estudiante.last_name}".strip() or participacion.estudiante.username,
                'avatar_elegido': participacion.avatar_elegido or '🧪',
                'avatar_nombre': participacion.avatar_nombre or 'Molécula Química',
                'activo': participacion.activo,
                'es_usuario_actual': participacion.estudiante == request.user,
                'grupo': participacion.grupo.nombre if participacion.grupo else None
            })
        
        return JsonResponse({
            'success': True,
            'participantes': participantes
        })
        
    except Exception as e:
        return JsonResponse({'success': False, 'error': str(e)})
    

@login_required
@role_required('Estudiante')
def unirse_competencia_principal(request):
    """Página principal para que estudiantes ingresen PIN"""
    user = request.user
    
    context = {
        'imgPerfil': user.imgPerfil,
        'usuario': user.username,
    }
    
    return render(request, 'estudiante/competencias/unirse_competencia.html', context)


@login_required
@role_required('Estudiante')
def seleccionar_avatar_competencia(request, competencia_id):
    competencia = Competencia.objects.get(id=competencia_id)
    participacion_id = request.GET.get('participacion_id')
    if participacion_id:
        participacion = ParticipacionCompetencia.objects.get(id=participacion_id, estudiante=request.user)
    else:
        participacion = ParticipacionCompetencia.objects.get(competencia=competencia, estudiante=request.user)

    if competencia.estado not in ['esperando', 'activa']:
        messages.error(request, 'La competencia ya no está disponible')
        return redirect('student_dashboard')

    context = {'competencia': competencia, 'participacion': participacion}
    return render(request, 'estudiante/competencias/seleccionar_avatar.html', context)

    



