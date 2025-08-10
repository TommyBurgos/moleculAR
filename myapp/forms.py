from django import forms
from user.models import Curso, Recurso, Rol

from django.contrib.auth import get_user_model

User = get_user_model()

class EditarPerfilForm(forms.ModelForm):
    class Meta:
        model = User
        fields = ['first_name', 'last_name', 'email', 'imgPerfil']
        widgets = {
            'first_name': forms.TextInput(attrs={'class': 'form-control', 'placeholder': 'Nombre'}),
            'last_name': forms.TextInput(attrs={'class': 'form-control', 'placeholder': 'Apellido'}),
            'email': forms.EmailInput(attrs={'class': 'form-control', 'placeholder': 'Correo electrónico'}),
        }

class CursoForm(forms.ModelForm):
    class Meta:
        model = Curso
        fields = ['titulo', 'descripcion', 'imagen_portada', 'categoria']
        widgets = {
            'titulo': forms.TextInput(attrs={'class': 'form-control'}),
            'descripcion': forms.Textarea(attrs={'class': 'form-control', 'rows': 4}),
            'categoria': forms.Select(attrs={'class': 'form-select'}),
        }

class RecursoForm(forms.ModelForm):
    class Meta:
        model = Recurso
        fields = ['titulo', 'descripcion', 'tipo', 'contenido_texto', 'video_url']
        widgets = {
            'contenido_texto': forms.Textarea(attrs={'rows': 4}),
            'video_url': forms.HiddenInput(),  # Se mostrará con JS al subir el video
        }

    def clean(self):
        cleaned_data = super().clean()
        tipo = cleaned_data.get('tipo')
        contenido = cleaned_data.get('contenido_texto')
        video = cleaned_data.get('video_url')

        if tipo and tipo.nombre.lower() == 'texto' and not contenido:
            raise forms.ValidationError("Debe ingresar el contenido del texto.")
        if tipo and tipo.nombre.lower() == 'video' and not video:
            raise forms.ValidationError("Debe subir un video.")

class UnirseCursoForm(forms.Form):
    codigo_acceso = forms.CharField(label="Código del curso", max_length=10)


class AdminCrearUsuarioForm(forms.Form):
    username = forms.CharField(
        max_length=150,
        widget=forms.TextInput(attrs={'class': 'form-control', 'placeholder': 'Nombre de usuario'})
    )
    first_name = forms.CharField(
        max_length=150, required=False,
        widget=forms.TextInput(attrs={'class': 'form-control', 'placeholder': 'Nombre'})
    )
    last_name = forms.CharField(
        max_length=150, required=False,
        widget=forms.TextInput(attrs={'class': 'form-control', 'placeholder': 'Apellido'})
    )
    email = forms.EmailField(
        widget=forms.EmailInput(attrs={'class': 'form-control', 'placeholder': 'Correo electrónico'})
    )
    rol = forms.ModelChoiceField(
        queryset=Rol.objects.all(),
        widget=forms.Select(attrs={'class': 'form-select'})
    )
    imgPerfil = forms.ImageField(required=False)
    password1 = forms.CharField(
        widget=forms.PasswordInput(attrs={'class': 'form-control', 'placeholder': 'Contraseña'})
    )
    password2 = forms.CharField(
        widget=forms.PasswordInput(attrs={'class': 'form-control', 'placeholder': 'Confirmar contraseña'})
    )
    is_active = forms.BooleanField(
        required=False,
        initial=True,
        widget=forms.CheckboxInput(attrs={'class': 'form-check-input'})
    )

    def clean_username(self):
        username = self.cleaned_data['username']
        if User.objects.filter(username=username).exists():
            raise forms.ValidationError("Este nombre de usuario ya existe.")
        return username

    def clean_email(self):
        email = self.cleaned_data['email']
        if User.objects.filter(email=email).exists():
            raise forms.ValidationError("Este correo ya está registrado.")
        return email

    def clean(self):
        cleaned = super().clean()
        p1 = cleaned.get('password1')
        p2 = cleaned.get('password2')
        if p1 and p2 and p1 != p2:
            self.add_error('password2', "Las contraseñas no coinciden.")
        return cleaned
