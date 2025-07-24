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


from user.models import User,Rol, Curso, Seccion, TipoRecurso, Recurso, Cuestionario, Practica, Modelo
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
