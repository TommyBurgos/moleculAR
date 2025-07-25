from django.shortcuts import render, redirect, get_object_or_404
from django.http import HttpResponse
from django.shortcuts import render
from django.http import JsonResponse
from django.views.decorators.csrf import csrf_exempt
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import SDWriter
from django.middleware.csrf import rotate_token
from django.contrib.auth.decorators import login_required
from django.contrib import messages
import tempfile
import json


from user.models import User,Rol, Curso, Seccion, TipoRecurso, Recurso, Cuestionario, Practica, PreguntaCuestionario, OpcionPregunta, TipoPregunta
from django.contrib.auth import login, logout, authenticate
import re
from .permissions import role_required

@csrf_exempt
def modificar_molecula(request):
    if request.method == 'POST':
        data = json.loads(request.body)
        smiles = data.get("smiles")
        atomo_idx = int(data.get("atomo_idx"))
        grupo = data.get("grupo")

        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        edmol = Chem.RWMol(mol)

        try:
            # Remover un hidrógeno unido al átomo objetivo
            h_idx = None
            for vecino in edmol.GetAtomWithIdx(atomo_idx).GetNeighbors():
                if vecino.GetAtomicNum() == 1:
                    h_idx = vecino.GetIdx()
                    break

            if h_idx is not None:
                edmol.RemoveAtom(h_idx)

            if grupo == "CH3":
                c = Chem.Atom(6)
                h1 = Chem.Atom(1)
                h2 = Chem.Atom(1)
                h3 = Chem.Atom(1)
                c_idx = edmol.AddAtom(c)
                edmol.AddBond(atomo_idx, c_idx, Chem.rdchem.BondType.SINGLE)
                edmol.AddBond(c_idx, edmol.AddAtom(h1), Chem.rdchem.BondType.SINGLE)
                edmol.AddBond(c_idx, edmol.AddAtom(h2), Chem.rdchem.BondType.SINGLE)
                edmol.AddBond(c_idx, edmol.AddAtom(h3), Chem.rdchem.BondType.SINGLE)

            mol_final = edmol.GetMol()
            Chem.SanitizeMol(mol_final)
            mol_final.UpdatePropertyCache(strict=False)
            AllChem.EmbedMolecule(mol_final)
            AllChem.UFFOptimizeMolecule(mol_final)

            # Guardar SDF temporal
            with tempfile.NamedTemporaryFile(suffix=".sdf", delete=False) as temp:
                writer = SDWriter(temp.name)
                writer.write(mol_final)
                writer.close()
                with open(temp.name, "r") as f:
                    sdf_data = f.read()

            return JsonResponse({"sdf": sdf_data})
        except Exception as e:
            return JsonResponse({"error": str(e)}, status=400)

    # ✅ Manejo para GET
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
                # Si quieres asignar rol por defecto
                rol_estudiante = Rol.objects.get(nombre='Estudiante')
                print(rol_estudiante)
                user.rol = rol_estudiante
                print(user.rol)                                
                user.save()
                print('***FIN***')
                login(request, user)
                print("Usuario creado y logueado exitosamente.")
                
                return redirect('inicio')  # O la vista que tengas definida
            except Exception as e:
                print(f"❌ Error al crear el usuario: {e}")
                return HttpResponse('Ocurrió un error al crear el usuario.')
        else:
            return HttpResponse('Las contraseñas no coinciden.')
    return HttpResponse('Método no permitido')


def vistaLogin(request):
    return render(request,'sign-in.html')

@role_required('Admin')
def inicioAdmin(request):
    user = request.user    
    imgPerfil=user.imgPerfil  
    context = {                  
        'imgPerfil': imgPerfil,        
        'usuario':user.username,        
    }   
    return render(request, 'usAdmin/admin-dashboard.html', context)

@role_required('Admin')
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

@role_required('Admin')
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

def ver_practica(request, recurso_id):
    recurso = get_object_or_404(Recurso, id=recurso_id)
    practica = get_object_or_404(Practica, recurso=recurso)
    print(f"el recurso es {recurso}")
    print(f"la practica es {practica}")

    modelo = practica.modelo_objetivo
    print(f"el modelo es {modelo}")
    print(f"el smile es {modelo.smiles}")
    print(modelo.smiles if modelo else "")
    smiles = modelo.smiles if modelo else ""
    
    return render(request, 'usAdmin/ver_practica.html', {
        'practica': practica,
        'recurso': recurso,
        'modelo':modelo,
        'smiles': smiles
    })


def detalle_recurso(request, recurso_id):
    recurso = get_object_or_404(Recurso, id=recurso_id)
    curso = recurso.seccion.curso 
    recursos_prefetch = Prefetch(
        'recursos',
        queryset=Recurso.objects.order_by('orden'),
        to_attr='recursos_ordenados'
    )    
    secciones = Seccion.objects.filter(curso=curso).prefetch_related(recursos_prefetch).order_by('orden')
    return render(request, 'usAdmin/detalleRecurso.html', {
        'recurso': recurso,
        'secciones': secciones,
        'curso': curso})


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
        response = requests.get(url, timeout=10)
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
def editar_practica(request, recurso_id):
    recurso = get_object_or_404(Recurso, id=recurso_id)
    practica = Practica.objects.filter(recurso=recurso).first()

    if request.method == 'POST':
        fuente = request.POST.get('fuente_molecula')
        titulo_objetivo = request.POST.get('titulo_objetivo')
        instrucciones = request.POST.get('instrucciones')
        smiles = None
        modelo = None

        if fuente == 'jsme':
            smiles = request.POST.get('mol_input')
            if smiles and smiles.strip():
                modelo = Modelo.objects.create(
                    titulo=titulo_objetivo,
                    descripcion=instrucciones,
                    tipo='2D',
                    smiles=smiles.strip(),
                    creado_por=request.user,
                    visible_biblioteca=False
                )
            else:
                return render(request, 'editar_practica.html', {
                    'recurso': recurso,
                    'error': 'No se pudo obtener el SMILES desde el editor JSME.'
                })

        elif fuente == 'pubchem':
            nombre_molecula = request.POST.get('nombre_molecula')
            print("Nombre recibido desde formulario:", nombre_molecula)
            smiles = obtener_smiles_pubchem(nombre_molecula)
            print("SMILES recibido desde PubChem:", smiles)
            if smiles and smiles.strip():
                modelo = Modelo.objects.create(
                    titulo=nombre_molecula,
                    descripcion='Molécula importada desde PubChem',
                    tipo='2D',
                    smiles=smiles.strip(),
                    creado_por=request.user,
                    visible_biblioteca=False
                )
            else:
                return render(request, 'usAdmin/editor_dibujo.html', {
                    'recurso': recurso,
                    'error': f"No se pudo encontrar un SMILES válido para '{nombre_molecula}'. Intenta con otro nombre."
                })

        # Crear la práctica solo si se creó un modelo
        if modelo:
            practica, creada = Practica.objects.get_or_create(
            recurso=recurso,
            defaults={
                'titulo_objetivo': titulo_objetivo,
                'instrucciones': instrucciones,
                'modelo_objetivo': modelo
            }
        )

        if not creada:
            practica.titulo_objetivo = titulo_objetivo
            practica.instrucciones = instrucciones
            practica.modelo_objetivo = modelo
            practica.save()

            return redirect('detalle_curso', curso_id=recurso.seccion.curso.id)
        else:
            return render(request, 'usAdmin/editor_dibujo.html', {
                'recurso': recurso,
                'error': 'No se pudo crear la práctica porque no se generó ningún modelo.'
            })

    return render(request, 'usAdmin/editor_dibujo.html', {'recurso': recurso})


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
        seccion_id = request.POST.get('seccion_id')
        print("La sección id es" + seccion_id)
        tipo = TipoRecurso.objects.get(id=tipo_id)
        seccion = Seccion.objects.get(id=seccion_id)

        recurso = Recurso.objects.create(
            titulo=titulo,
            descripcion=descripcion,
            tipo=tipo,
            imagen=imagen,
            seccion=seccion,
            visible_biblioteca=visible
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


        elif 'cuestionario' in nombre_tipo:
            print("Dentro del if cuestionario")
            instrucciones = request.POST.get('nombre_cuestionario')
            tiempo_limite = request.POST.get('descripcion_cuestionario')
            Cuestionario.objects.create(
                recurso=recurso,
                instrucciones=instrucciones,
                tiempo_limite=tiempo_limite
            )        
        elif tipo.nombre.lower() == 'practica':
            # Por ahora NO creamos la práctica. Solo dejamos el recurso creado.
            # Luego redirigiremos a una vista donde se edita ese recurso.
            print("Se creo el recurso practica")
            return redirect('editar_practica', recurso.id)                





import boto3
import uuid
from django.conf import settings

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

        # ⚠️ NO incluir campos ACL ni condiciones ACL
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

@role_required('Admin')
def vistaCrearCurso(request):
    user = request.user    
    imgPerfil=user.imgPerfil 
    if request.method == 'POST':
        form = CursoForm(request.POST, request.FILES)
        if form.is_valid():
            curso = form.save(commit=False)
            curso.profesor = request.user  # Asignar el usuario actual como profesor
            curso.save()
            return redirect('detalleCursos-adm', curso.id) 
    else:
        form = CursoForm() 
    context = {                  
        'imgPerfil': imgPerfil,        
        'usuario':user.username,
        'form': form        
    }   
    return render(request,'usAdmin/crearCurso.html',context)



@role_required('Estudiante')
def inicioEstudiante(request):
    user = request.user    
    return render(request, 'estudiante/student-dashboard.html')    

def inicioDocente(request):
    return render(request, 'docente/instructor-dashboard.html')



def custom_login(request):
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
            if rol == 1:
                print("Ingrese al rol 1")
                return redirect('dashboard-adm')
            elif rol == 2:
                return redirect('student_dashboard')
            elif rol == 3:
                return redirect('teacher_dashboard')
        else:
            # Manejar error de autenticación
            return render(request, 'generales/accesoDenegado.html', {'error': 'Credenciales inválidas'})
    return redirect('inicio')

def signout(request):
    logout(request)
    rotate_token(request)  # Gira el token CSRF para la nueva sesión
    return redirect('inicio')

from rdkit import Chem
from rdkit.Chem import AllChem, SDWriter
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
            if tipo == 'smiles':
                mol = Chem.MolFromSmiles(entrada)
                mol = Chem.AddHs(mol)
            elif tipo == 'mol':
                mol = Chem.MolFromMolBlock(entrada)
            else:
                return JsonResponse({'error': 'Formato no válido'}, status=400)

            AllChem.EmbedMolecule(mol)
            AllChem.UFFOptimizeMolecule(mol)

            with tempfile.NamedTemporaryFile(suffix=".sdf", delete=False) as temp:
                writer = SDWriter(temp.name)
                writer.write(mol)
                writer.close()
                with open(temp.name, "r") as f:
                    sdf_data = f.read()

            return render(request, 'visor_3d.html', {'sdf_data': sdf_data})

        except Exception as e:
            return JsonResponse({'error': str(e)}, status=400)

    return JsonResponse({'message': 'Usa POST con SMILES o MOL'})


#Paula

from django.views.decorators.http import require_http_methods
from django.utils import timezone
from django.db import transaction
from django.db.models import Sum, Max
from datetime import datetime

#@login_required
def instructorQuiz(request, seccion_id=None):
    """Vista principal para crear cuestionarios - MEJORADA"""
    seccion = None
    if seccion_id:
        seccion = get_object_or_404(Seccion, id=seccion_id)
        # Verificar que el usuario sea el profesor del curso
        if seccion.curso.profesor != request.user:
            return redirect('dashboard')
    
    # CAMBIO: No crear cuestionario automáticamente en POST
    # Solo mostrar el formulario, el cuestionario se crea al hacer clic en "Guardar"
    
    # NUEVO: Obtener tipos de pregunta desde la base de datos
    tipos_pregunta = TipoPregunta.objects.all()
    
    # Si no hay tipos de pregunta, crear los básicos
    if not tipos_pregunta.exists():
        tipos_pregunta = crear_tipos_pregunta_basicos()
    
    context = {
        'cuestionario': None,
        'preguntas': [],
        'seccion': seccion,
        'seccion_id': seccion_id,
        'tipos_pregunta': tipos_pregunta,  # NUEVO: Pasar tipos de pregunta al template
    }
    return render(request, 'docente/crear_cuestionario.html', context)

# NUEVA FUNCIÓN AUXILIAR: Crear tipos de pregunta básicos
def crear_tipos_pregunta_basicos():
    """Crear tipos de pregunta básicos si no existen"""
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
# 2. FUNCIÓN EDITAR CUESTIONARIO (consistencia con tipos dinámicos)
# =====================================================

#@login_required
def editarCuestionario(request, cuestionario_id):
    """Vista para editar un cuestionario existente"""
    cuestionario = get_object_or_404(Cuestionario, id=cuestionario_id)
    
    # Verificar permisos - solo si hay sección
    if cuestionario.recurso.seccion and cuestionario.recurso.seccion.curso.profesor != request.user:
        return redirect('dashboard')
    
    preguntas = cuestionario.preguntas.all().order_by('orden')
    tipos_pregunta = TipoPregunta.objects.all()  # CONSISTENCIA: Obtener tipos dinámicamente
    
    # Calcular puntaje total y verificar calificación automática
    puntaje_calculado = preguntas.aggregate(total=Sum('puntaje'))['total'] or 0
    tiene_preguntas_manuales = preguntas.filter(
        tipo__nombre__in=['respuesta_abierta', 'simulador_2d', 'simulador_3d']
    ).exists()
    
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
# 3. FUNCIÓN ACTUALIZAR CONFIGURACIÓN MEJORADA
# =====================================================

@csrf_exempt
@require_http_methods(["POST"])
#@login_required
def actualizarConfiguracion(request):
    """Actualizar configuración general del cuestionario - MEJORADA"""
    try:
        data = json.loads(request.body)
        cuestionario_id = data.get('cuestionario_id')
        configuracion = data.get('configuracion')  # NUEVO: para el modal de configuración
        
        cuestionario = get_object_or_404(Cuestionario, id=cuestionario_id)
        
        # Verificar permisos (solo si hay sección)
        if cuestionario.recurso.seccion and cuestionario.recurso.seccion.curso.profesor != request.user:
            return JsonResponse({'success': False, 'error': 'Sin permisos'}, status=403)
        
        # NUEVO: Manejar datos del modal de configuración
        if configuracion:
            # Configuración avanzada del modal
            if 'tiempo_limite' in configuracion:
                cuestionario.tiempo_limite = int(configuracion['tiempo_limite'])
            
            if 'puntos_totales' in configuracion:
                # El puntaje se calcula automáticamente, pero se puede permitir override
                pass
            
            # NUEVO: Campos adicionales si existen en tu modelo
            if hasattr(cuestionario, 'intentos_permitidos') and 'intentos_permitidos' in configuracion:
                intentos = configuracion['intentos_permitidos']
                cuestionario.intentos_permitidos = 999 if intentos == 'ilimitado' else int(intentos)
            
            if hasattr(cuestionario, 'mostrar_resultados') and 'mostrar_resultados' in configuracion:
                cuestionario.mostrar_resultados = bool(configuracion['mostrar_resultados'])
            
            if hasattr(cuestionario, 'orden_aleatorio') and 'mezclar_preguntas' in configuracion:
                cuestionario.orden_aleatorio = bool(configuracion['mezclar_preguntas'])
            
            if hasattr(cuestionario, 'tipo_calificacion') and 'tipo_calificacion' in configuracion:
                cuestionario.tipo_calificacion = configuracion['tipo_calificacion']
            
            # NUEVO: Manejar fechas del modal
            if hasattr(cuestionario, 'fecha_apertura') and 'fecha_inicio' in configuracion and configuracion['fecha_inicio']:
                try:
                    cuestionario.fecha_apertura = datetime.fromisoformat(
                        configuracion['fecha_inicio'].replace('Z', '+00:00')
                    )
                except:
                    pass
            
            if hasattr(cuestionario, 'fecha_cierre') and 'fecha_fin' in configuracion and configuracion['fecha_fin']:
                try:
                    cuestionario.fecha_cierre = datetime.fromisoformat(
                        configuracion['fecha_fin'].replace('Z', '+00:00')
                    )
                except:
                    pass
        else:
            # Configuración básica (código original)
            if 'titulo' in data:
                cuestionario.recurso.titulo = data['titulo']
                cuestionario.recurso.save()
            
            if 'instrucciones' in data:
                cuestionario.instrucciones = data['instrucciones']
            
            if 'tiempo_limite' in data:
                cuestionario.tiempo_limite = int(data['tiempo_limite'])
            
            # Tu código original para fechas
            if hasattr(cuestionario, 'fecha_apertura') and 'fecha_apertura' in data and data['fecha_apertura']:
                try:
                    cuestionario.fecha_apertura = datetime.fromisoformat(
                        data['fecha_apertura'].replace('Z', '+00:00')
                    )
                except:
                    pass
            
            if hasattr(cuestionario, 'fecha_cierre') and 'fecha_cierre' in data and data['fecha_cierre']:
                try:
                    cuestionario.fecha_cierre = datetime.fromisoformat(
                        data['fecha_cierre'].replace('Z', '+00:00')
                    )
                except:
                    pass
            
            # Resto de tu código original...
            if hasattr(cuestionario, 'intentos_permitidos') and 'intentos_permitidos' in data:
                cuestionario.intentos_permitidos = int(data['intentos_permitidos'])
            
            if hasattr(cuestionario, 'mostrar_resultados') and 'mostrar_resultados' in data:
                cuestionario.mostrar_resultados = bool(data['mostrar_resultados'])
            
            if hasattr(cuestionario, 'orden_aleatorio') and 'orden_aleatorio' in data:
                cuestionario.orden_aleatorio = bool(data['orden_aleatorio'])
        
        # Calcular puntaje total automáticamente (código original)
        puntaje_total = cuestionario.preguntas.aggregate(
            total=Sum('puntaje')
        )['total'] or 0
        
        if hasattr(cuestionario, 'puntaje_total'):
            cuestionario.puntaje_total = puntaje_total
        
        # Determinar si requiere calificación manual (código original)
        tiene_preguntas_manuales = cuestionario.preguntas.filter(
            tipo__nombre__in=['respuesta_abierta', 'simulador_2d', 'simulador_3d']
        ).exists()
        
        if hasattr(cuestionario, 'calificacion_automatica'):
            cuestionario.calificacion_automatica = not tiene_preguntas_manuales
        
        cuestionario.save()
        
        return JsonResponse({
            'success': True,
            'mensaje': 'Configuración actualizada exitosamente',
            'puntaje_total': float(puntaje_total),
            'calificacion_automatica': not tiene_preguntas_manuales
        })
        
    except Exception as e:
        return JsonResponse({
            'success': False,
            'error': str(e)
        }, status=400)

# =====================================================
# 4. FUNCIÓN AGREGAR PREGUNTA (mantén tu función original)
# =====================================================

@csrf_exempt
@require_http_methods(["POST"])
#@login_required
def agregarPregunta(request):
    """Agregar una nueva pregunta al cuestionario"""
    try:
        data = json.loads(request.body)
        cuestionario_id = data.get('cuestionario_id')
        tipo_nombre = data.get('tipo')
        
        cuestionario = get_object_or_404(Cuestionario, id=cuestionario_id)
        
        # Verificar permisos (solo si hay sección)
        if cuestionario.recurso.seccion and cuestionario.recurso.seccion.curso.profesor != request.user:
            return JsonResponse({'success': False, 'error': 'Sin permisos'}, status=403)
        
        # Obtener el tipo de pregunta
        tipo_pregunta = get_object_or_404(TipoPregunta, nombre=tipo_nombre)
        
        # Obtener el próximo número de orden
        ultima_pregunta = cuestionario.preguntas.order_by('-orden').first()
        orden = ultima_pregunta.orden + 1 if ultima_pregunta else 1
        
        pregunta = PreguntaCuestionario.objects.create(
            cuestionario=cuestionario,
            tipo=tipo_pregunta,
            enunciado=f"Nueva pregunta {orden}",
            orden=orden,
            puntaje=1
        )
        
        # Crear opciones por defecto para preguntas de opción única/múltiple
        if tipo_nombre in ['opcion_unica', 'opcion_multiple', 'seleccion_multiple', 'falso_verdadero']:
            if tipo_nombre == 'falso_verdadero':
                opciones_default = [
                    ('Verdadero', True),
                    ('Falso', False)
                ]
            else:
                opciones_default = [
                    ('Opción 1', True),
                    ('Opción 2', False)
                ]
            
            for texto, es_correcta in opciones_default:
                OpcionPregunta.objects.create(
                    pregunta=pregunta,
                    texto=texto,
                    es_correcta=es_correcta
                )
        
        # Recalcular puntaje total del cuestionario
        puntaje_total = cuestionario.preguntas.aggregate(
            total=Sum('puntaje')
        )['total'] or 0
        
        if hasattr(cuestionario, 'puntaje_total'):
            cuestionario.puntaje_total = puntaje_total
        
        # Actualizar calificación automática
        tiene_preguntas_manuales = cuestionario.preguntas.filter(
            tipo__nombre__in=['respuesta_abierta', 'simulador_2d', 'simulador_3d']
        ).exists()
        
        if hasattr(cuestionario, 'calificacion_automatica'):
            cuestionario.calificacion_automatica = not tiene_preguntas_manuales
            
        cuestionario.save()
        
        return JsonResponse({
            'success': True,
            'pregunta_id': pregunta.id,
            'mensaje': 'Pregunta agregada exitosamente',
            'puntaje_total': float(puntaje_total),
            'calificacion_automatica': not tiene_preguntas_manuales
        })
        
    except Exception as e:
        return JsonResponse({
            'success': False,
            'error': str(e)
        }, status=400)

# =====================================================
# 5. FUNCIÓN AGREGAR OPCIÓN MEJORADA
# =====================================================

@csrf_exempt
@require_http_methods(["POST"])
#@login_required
def agregarOpcion(request):
    """Agregar una nueva opción a una pregunta - MEJORADA"""
    try:
        data = json.loads(request.body)
        pregunta_id = data.get('pregunta_id')
        
        pregunta = get_object_or_404(PreguntaCuestionario, id=pregunta_id)
        
        # Verificar permisos (solo si hay sección)
        if pregunta.cuestionario.recurso.seccion and pregunta.cuestionario.recurso.seccion.curso.profesor != request.user:
            return JsonResponse({'success': False, 'error': 'Sin permisos'}, status=403)
        
        # NUEVO: Verificar que el tipo de pregunta permita agregar opciones
        if pregunta.tipo.nombre not in ['opcion_unica', 'opcion_multiple']:
            return JsonResponse({
                'success': False,
                'error': 'No se pueden agregar opciones a este tipo de pregunta'
            })
        
        # Obtener el próximo número de opción
        opciones_count = pregunta.opciones.count()
        
        opcion = OpcionPregunta.objects.create(
            pregunta=pregunta,
            texto=f"Nueva opción {opciones_count + 1}",
            es_correcta=False
        )
        
        return JsonResponse({
            'success': True,
            'opcion_id': opcion.id,
            'mensaje': 'Opción agregada exitosamente',
            'mantener_acordeon': True  # NUEVO: Flag para mantener acordeón abierto
        })
        
    except Exception as e:
        return JsonResponse({
            'success': False,
            'error': str(e)
        }, status=400)

# =====================================================
# 6. FUNCIÓN ACTUALIZAR OPCIÓN (mantén tu función original)
# =====================================================

@csrf_exempt
@require_http_methods(["POST"])
#@login_required
def actualizarOpcion(request):
    """Actualizar una opción existente"""
    try:
        data = json.loads(request.body)
        opcion_id = data.get('opcion_id')
        texto = data.get('texto')
        es_correcta = data.get('es_correcta', False)
        
        opcion = get_object_or_404(OpcionPregunta, id=opcion_id)
        
        # Verificar permisos (solo si hay sección)
        if opcion.pregunta.cuestionario.recurso.seccion and opcion.pregunta.cuestionario.recurso.seccion.curso.profesor != request.user:
            return JsonResponse({'success': False, 'error': 'Sin permisos'}, status=403)
        
        opcion.texto = texto
        opcion.es_correcta = es_correcta
        opcion.save()
        
        return JsonResponse({
            'success': True,
            'mensaje': 'Opción actualizada exitosamente'
        })
        
    except Exception as e:
        return JsonResponse({
            'success': False,
            'error': str(e)
        }, status=400)

# =====================================================
# 7. FUNCIÓN ACTUALIZAR PREGUNTA MEJORADA
# =====================================================

@csrf_exempt
@require_http_methods(["POST"])
#@login_required
def actualizarPregunta(request):
    """Actualizar el enunciado de una pregunta - MEJORADA"""
    try:
        data = json.loads(request.body)
        pregunta_id = data.get('pregunta_id')
        enunciado = data.get('enunciado')
        puntaje = data.get('puntaje', 1)
        calificacion_manual = data.get('calificacion_manual', False)  # NUEVO
        
        pregunta = get_object_or_404(PreguntaCuestionario, id=pregunta_id)
        
        # Verificar permisos (solo si hay sección)
        if pregunta.cuestionario.recurso.seccion and pregunta.cuestionario.recurso.seccion.curso.profesor != request.user:
            return JsonResponse({'success': False, 'error': 'Sin permisos'}, status=403)
        
        pregunta.enunciado = enunciado
        pregunta.puntaje = int(puntaje)
        
        # NUEVO: Manejar calificación manual si el campo existe
        if hasattr(pregunta, 'calificacion_manual'):
            pregunta.calificacion_manual = bool(calificacion_manual)
        
        pregunta.save()
        
        # Recalcular puntaje total del cuestionario
        cuestionario = pregunta.cuestionario
        puntaje_total = cuestionario.preguntas.aggregate(
            total=Sum('puntaje')
        )['total'] or 0
        
        if hasattr(cuestionario, 'puntaje_total'):
            cuestionario.puntaje_total = puntaje_total
            cuestionario.save()
        
        return JsonResponse({
            'success': True,
            'mensaje': 'Pregunta actualizada exitosamente',
            'puntaje_total': float(puntaje_total)
        })
        
    except Exception as e:
        return JsonResponse({
            'success': False,
            'error': str(e)
        }, status=400)

# =====================================================
# 8. FUNCIÓN ELIMINAR OPCIÓN (mantén tu función original)
# =====================================================

@csrf_exempt
@require_http_methods(["DELETE"])
#@login_required
def eliminarOpcion(request, opcion_id):
    """Eliminar una opción"""
    try:
        opcion = get_object_or_404(OpcionPregunta, id=opcion_id)
        
        # Verificar permisos (solo si hay sección)
        if opcion.pregunta.cuestionario.recurso.seccion and opcion.pregunta.cuestionario.recurso.seccion.curso.profesor != request.user:
            return JsonResponse({'success': False, 'error': 'Sin permisos'}, status=403)
        
        # No permitir eliminar si solo quedan 2 opciones
        if opcion.pregunta.opciones.count() <= 2:
            return JsonResponse({
                'success': False,
                'error': 'No se puede eliminar. Mínimo 2 opciones requeridas.'
            }, status=400)
        
        opcion.delete()
        
        return JsonResponse({
            'success': True,
            'mensaje': 'Opción eliminada exitosamente'
        })
        
    except Exception as e:
        return JsonResponse({
            'success': False,
            'error': str(e)
        }, status=400)

# =====================================================
# 9. FUNCIÓN ELIMINAR PREGUNTA (mantén tu función original)
# =====================================================

@csrf_exempt
@require_http_methods(["DELETE"])
#@login_required
def eliminarPregunta(request, pregunta_id):
    """Eliminar una pregunta completa"""
    try:
        pregunta = get_object_or_404(PreguntaCuestionario, id=pregunta_id)
        
        # Verificar permisos (solo si hay sección)
        if pregunta.cuestionario.recurso.seccion and pregunta.cuestionario.recurso.seccion.curso.profesor != request.user:
            return JsonResponse({'success': False, 'error': 'Sin permisos'}, status=403)
        
        cuestionario = pregunta.cuestionario
        pregunta.delete()
        
        # Recalcular puntaje total del cuestionario
        puntaje_total = cuestionario.preguntas.aggregate(
            total=Sum('puntaje')
        )['total'] or 0
        
        if hasattr(cuestionario, 'puntaje_total'):
            cuestionario.puntaje_total = puntaje_total
        
        # Actualizar calificación automática
        tiene_preguntas_manuales = cuestionario.preguntas.filter(
            tipo__nombre__in=['respuesta_abierta', 'simulador_2d', 'simulador_3d']
        ).exists()
        
        if hasattr(cuestionario, 'calificacion_automatica'):
            cuestionario.calificacion_automatica = not tiene_preguntas_manuales
            
        cuestionario.save()
        
        return JsonResponse({
            'success': True,
            'mensaje': 'Pregunta eliminada exitosamente',
            'puntaje_total': float(puntaje_total),
            'calificacion_automatica': not tiene_preguntas_manuales
        })
        
    except Exception as e:
        return JsonResponse({
            'success': False,
            'error': str(e)
        }, status=400)

# =====================================================
# 10. FUNCIÓN CREAR TIPOS PREGUNTA (mantén tu función original)
# =====================================================

@csrf_exempt
@require_http_methods(["POST"])
#@login_required
def crearTiposPregunta(request):
    """Crear tipos de pregunta automáticamente"""
    try:
        # Crear tipos de pregunta
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
        
        # Crear tipo de recurso cuestionario
        tipo_recurso, created = TipoRecurso.objects.get_or_create(
            nombre='cuestionario',
            defaults={'descripcion': 'Recurso tipo cuestionario para evaluaciones'}
        )
        
        return JsonResponse({
            'success': True,
            'mensaje': f'Tipos de pregunta listos. {creados} nuevos creados.',
        })
        
    except Exception as e:
        return JsonResponse({
            'success': False,
            'error': str(e)
        }, status=400)

# =====================================================
# 11. NUEVA FUNCIÓN: CAMBIAR TIPO DE PREGUNTA
# =====================================================

@csrf_exempt
@require_http_methods(["POST"])
#@login_required
def cambiarTipoPregunta(request):
    """NUEVA: Cambiar el tipo de una pregunta existente"""
    try:
        data = json.loads(request.body)
        pregunta_id = data.get('pregunta_id')
        nuevo_tipo = data.get('nuevo_tipo')
        
        pregunta = get_object_or_404(PreguntaCuestionario, id=pregunta_id)
        
        # Verificar permisos
        if pregunta.cuestionario.recurso.seccion and pregunta.cuestionario.recurso.seccion.curso.profesor != request.user:
            return JsonResponse({'success': False, 'error': 'Sin permisos'}, status=403)
        
        with transaction.atomic():
            # Eliminar opciones existentes
            pregunta.opciones.all().delete()
            
            # Cambiar tipo
            tipo_pregunta = get_object_or_404(TipoPregunta, nombre=nuevo_tipo)
            pregunta.tipo = tipo_pregunta
            
            # Resetear calificación manual si el campo existe
            if hasattr(pregunta, 'calificacion_manual'):
                pregunta.calificacion_manual = False
            
            pregunta.save()
            
            # Crear nuevas opciones según el tipo
            if nuevo_tipo in ['opcion_unica', 'opcion_multiple']:
                opciones_default = [
                    ('Opción 1', True),
                    ('Opción 2', False),
                    ('Opción 3', False),
                    ('Opción 4', False)
                ]
                for texto, es_correcta in opciones_default:
                    OpcionPregunta.objects.create(
                        pregunta=pregunta,
                        texto=texto,
                        es_correcta=es_correcta
                    )
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
        
        return JsonResponse({
            'success': True,
            'mensaje': 'Tipo de pregunta cambiado exitosamente'
        })
        
    except Exception as e:
        return JsonResponse({
            'success': False,
            'error': str(e)
        }, status=400)

# =====================================================
# 12. NUEVA FUNCIÓN: GUARDAR CUESTIONARIO COMPLETO CON PREGUNTAS Y CONFIGURACIÓN
# =====================================================

@csrf_exempt
@require_http_methods(["POST"])
#@login_required
def guardarCuestionarioCompleto(request):
    """NUEVA: Guardar cuestionario completo con preguntas temporales y configuración"""
    try:
        # Obtener datos básicos
        seccion_id = request.POST.get('seccion_id')
        titulo = request.POST.get('titulo', '').strip()
        descripcion = request.POST.get('descripcion', '').strip()
        
        # NUEVO: Obtener preguntas temporales y configuración
        preguntas_temporales_json = request.POST.get('preguntas_temporales', '[]')
        configuracion_temporal_json = request.POST.get('configuracion_temporal', '{}')
        
        try:
            preguntas_temporales = json.loads(preguntas_temporales_json)
        except:
            preguntas_temporales = []
            
        try:
            configuracion_temporal = json.loads(configuracion_temporal_json)
        except:
            configuracion_temporal = {}
        
        if not titulo:
            return JsonResponse({
                'success': False,
                'error': 'El título del cuestionario es obligatorio'
            })
        
        seccion = None
        if seccion_id:
            seccion = get_object_or_404(Seccion, id=seccion_id)
            # Verificar permisos
            if seccion.curso.profesor != request.user:
                return JsonResponse({'success': False, 'error': 'Sin permisos'}, status=403)
        
        with transaction.atomic():
            # Obtener o crear el tipo de recurso "cuestionario"
            tipo_cuestionario, created = TipoRecurso.objects.get_or_create(
                nombre='cuestionario',
                defaults={'descripcion': 'Recurso tipo cuestionario'}
            )
            
            # Crear el recurso
            recurso = Recurso.objects.create(
                seccion=seccion,
                tipo=tipo_cuestionario,
                titulo=titulo,
                descripcion=descripcion,
                orden=seccion.recursos.count() + 1 if seccion else 1
            )
            
            # Crear el cuestionario vinculado al recurso CON configuración temporal
            cuestionario_data = {
                'recurso': recurso,
                'instrucciones': descripcion,
                'tiempo_limite': configuracion_temporal.get('tiempo_limite', 30),
            }
            
            # Agregar campos de configuración si existen en el modelo
            if hasattr(Cuestionario, 'tipo_calificacion'):
                cuestionario_data['tipo_calificacion'] = configuracion_temporal.get('tipo_calificacion', 'automatica')
            if hasattr(Cuestionario, 'intentos_permitidos'):
                intentos = configuracion_temporal.get('intentos_permitidos', 1)
                cuestionario_data['intentos_permitidos'] = 999 if intentos == 'ilimitado' else int(intentos)
            if hasattr(Cuestionario, 'mostrar_resultados'):
                cuestionario_data['mostrar_resultados'] = configuracion_temporal.get('mostrar_resultados', True)
            if hasattr(Cuestionario, 'orden_aleatorio'):
                cuestionario_data['orden_aleatorio'] = configuracion_temporal.get('mezclar_preguntas', False)
            
            # Manejar fechas si existen en el modelo
            if hasattr(Cuestionario, 'fecha_apertura') and configuracion_temporal.get('fecha_inicio'):
                try:
                    cuestionario_data['fecha_apertura'] = datetime.fromisoformat(
                        configuracion_temporal['fecha_inicio'].replace('Z', '+00:00')
                    )
                except:
                    pass
                    
            if hasattr(Cuestionario, 'fecha_cierre') and configuracion_temporal.get('fecha_fin'):
                try:
                    cuestionario_data['fecha_cierre'] = datetime.fromisoformat(
                        configuracion_temporal['fecha_fin'].replace('Z', '+00:00')
                    )
                except:
                    pass
            
            # Crear el cuestionario con toda la configuración
            cuestionario = Cuestionario.objects.create(**cuestionario_data)
            
            # Manejar imagen si se subió
            if request.FILES.get('imagen') and hasattr(cuestionario, 'imagen'):
                cuestionario.imagen = request.FILES['imagen']
                cuestionario.save()
            
            # NUEVO: Crear preguntas temporales en la base de datos
            preguntas_creadas = 0
            for pregunta_temp in preguntas_temporales:
                # Obtener o crear el tipo de pregunta
                tipo_pregunta, created = TipoPregunta.objects.get_or_create(
                    nombre=pregunta_temp['tipo'],
                    defaults={'descripcion': f'Pregunta de tipo {pregunta_temp["tipo"]}'}
                )
                
                # Crear la pregunta
                pregunta_data = {
                    'cuestionario': cuestionario,
                    'tipo': tipo_pregunta,
                    'enunciado': pregunta_temp['enunciado'],
                    'orden': pregunta_temp['orden'],
                    'puntaje': int(pregunta_temp['puntaje'])
                }
                
                # Agregar calificación manual si el campo existe
                if hasattr(PreguntaCuestionario, 'calificacion_manual'):
                    pregunta_data['calificacion_manual'] = pregunta_temp.get('calificacion_manual', False)
                
                pregunta = PreguntaCuestionario.objects.create(**pregunta_data)
                
                # Crear las opciones si existen
                for opcion_temp in pregunta_temp.get('opciones', []):
                    OpcionPregunta.objects.create(
                        pregunta=pregunta,
                        texto=opcion_temp['texto'],
                        es_correcta=opcion_temp['es_correcta']
                    )
                
                preguntas_creadas += 1
            
            # Calcular puntaje total
            puntaje_total = cuestionario.preguntas.aggregate(
                total=Sum('puntaje')
            )['total'] or 0
            
            if hasattr(cuestionario, 'puntaje_total'):
                cuestionario.puntaje_total = puntaje_total
                cuestionario.save()
        
        # URL de redirección corregida con tu patrón de URLs
        redirect_url = f'/docente/editar_cuestionario/{cuestionario.id}/'
        
        return JsonResponse({
            'success': True,
            'cuestionario_id': cuestionario.id,
            'preguntas_creadas': preguntas_creadas,
            'puntaje_total': float(puntaje_total),
            'mensaje': f'Cuestionario creado exitosamente con {preguntas_creadas} preguntas',
            'redirect_url': redirect_url
        })
        
    except Exception as e:
        return JsonResponse({
            'success': False,
            'error': str(e)
        }, status=400)
