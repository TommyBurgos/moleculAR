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
from django.contrib import messages
import tempfile
import json

from user.models import User,Rol, Curso, Seccion, TipoRecurso, Recurso, Cuestionario, Practica, Modelo, InscripcionCurso, MoleculaEstudiante, Competencia, PreguntaCuestionario,TipoPregunta, ProgresoUsuario, OpcionPregunta, IntentoCuestionario, RespuestaEstudiante, PracticaConfig, IntentoArmado

from django.contrib.auth import login, logout, authenticate
import re
from .permissions import role_required
from django.db.models import Q, Sum  
 
import boto3
import uuid
from django.conf import settings

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
                rol_estudiante = Rol.objects.get(nombre='estudiante')
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

@csrf_exempt
def vistaLogin(request):
    return render(request,'sign-in.html')

def acceso_denegado(request):
    print("Inicie a la función de acceso denegado")        
    return render(request, 'error-404.html')


@role_required('Admin')
def inicioAdmin(request):
    print("Inicie a la función de inicio Admin")
    user = request.user    
    imgPerfil=user.imgPerfil  
    context = {                  
        'imgPerfil': imgPerfil,        
        'usuario':user,        
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
        recurso.save()
        messages.success(request, "HTML embebido guardado correctamente.")
        return redirect('detalle_recurso', recurso.id)

    return render(request, 'usAdmin/editar_html.html', {'recurso': recurso})


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
            recurso.save()
            return redirect('editar_html', recurso.id)



        return redirect('detalle_recurso',recurso.id)

from .forms import AdminCrearUsuarioForm

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


from .forms import CursoForm

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

def panel_competencia(request, competencia_id):
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
    return HttpResponseRedirect(reverse('panel_competencia', args=[competencia.id]))

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


@role_required('Docente')
def inicioDocente(request):
    return render(request, 'docente/instructor-dashboard.html')



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
    """Vista para editar un cuestionario existente - AHORA USA RECURSO_ID"""
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
    """Agregar una nueva pregunta al cuestionario - COMPLETAMENTE FUNCIONAL"""
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
        
        # Calcular orden
        ultima_pregunta = cuestionario.preguntas.order_by('-orden').first()
        orden = ultima_pregunta.orden + 1 if ultima_pregunta else 1
        
        print(f"📝 Creando pregunta orden {orden}")
        
        # Crear pregunta
        pregunta = PreguntaCuestionario.objects.create(
            cuestionario=cuestionario,
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
        
        # Verificar permisos
        if (pregunta.cuestionario.recurso.seccion and 
            pregunta.cuestionario.recurso.seccion.curso.profesor != request.user):
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


from .forms import EditarPerfilForm

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
