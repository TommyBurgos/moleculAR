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
