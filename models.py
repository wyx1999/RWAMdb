from django.db import models


class Site(models.Model):
    objects = models.Manager()
    id = models.CharField(max_length=100, primary_key=True, verbose_name="id")
    name = models.CharField(max_length=200, verbose_name="name")
    longitude = models.FloatField(null=True,blank=True, verbose_name="longitude")
    latitude = models.FloatField(null=True, blank=True, verbose_name="latitude")
    altitude = models.FloatField(null=True, blank=True, verbose_name="altitude")
    temperature = models.FloatField(null=True, blank=True, verbose_name="temperature")
    precipitation = models.FloatField(null=True, blank=True, verbose_name="precipitation")
    temperature_2022 = models.FloatField(null=True, blank=True, verbose_name="temperature_2022")
    precipitation_2022 = models.FloatField(null=True, blank=True, verbose_name="precipitation_2022")
    type = models.CharField(max_length=50, null=True, blank=True, verbose_name="type")
    location = models.CharField(max_length=200, verbose_name="location")
    introduction = models.TextField(null=True, blank=True, verbose_name="introduction")
    chinese = models.CharField(max_length=200, null=True, blank=True, verbose_name="chinese")
    old = models.CharField(max_length=200, null=True, blank=True, verbose_name="old")


class Sample(models.Model):
    objects = models.Manager()
    id = models.CharField(max_length=100, primary_key=True, verbose_name="id")
    site = models.ForeignKey(Site, on_delete=models.CASCADE, verbose_name="site", related_name='samples')

    TC = models.FloatField(null=True, verbose_name="TC")
    TS = models.FloatField(null=True, verbose_name="TS")
    TN = models.FloatField(null=True, verbose_name="TN")
    TP = models.FloatField(null=True, verbose_name="TP")
    SOC = models.FloatField(null=True, verbose_name="SOC")
    C_N = models.FloatField(null=True, verbose_name="C_N")
    IC = models.FloatField(null=True, verbose_name="IC")
    MBC = models.FloatField(null=True, verbose_name="MBC")
    pH = models.FloatField(null=True, verbose_name="pH")
    water = models.FloatField(null=True, verbose_name="water")
    NO3 = models.FloatField(null=True, verbose_name="NO3")
    NH4 = models.FloatField(null=True, verbose_name="NH4")

    date = models.CharField(max_length=100, null=True, verbose_name="date")
    person = models.CharField(max_length=200, null=True, verbose_name="person")
    old = models.CharField(max_length=200, null=True, blank=True, verbose_name="old")


class MAG(models.Model):
    objects = models.Manager()
    id = models.CharField(max_length=100, primary_key=True, verbose_name="id")
    sample = models.ForeignKey(Sample, on_delete=models.CASCADE, verbose_name="sample", related_name='MAGs')
    Kingdom = models.CharField(max_length=50, null=True)
    Phylum = models.CharField(max_length=50, null=True)
    Class = models.CharField(max_length=50, null=True)
    Order = models.CharField(max_length=50, null=True)
    Family = models.CharField(max_length=50, null=True)
    Genus = models.CharField(max_length=50, null=True)
    Species = models.CharField(max_length=50, null=True)
    taxonomy = models.CharField(max_length=200, null=True)
    relative_abundance = models.FloatField(null=True, verbose_name="relative_abundance")
    length = models.FloatField(null=True, verbose_name="length")
    RPKM = models.FloatField(null=True, verbose_name="RPKM")
    completeness = models.FloatField(null=True, verbose_name="Completeness")
    contamination = models.FloatField(null=True, verbose_name="Contamination")
    GC = models.FloatField(null=True, verbose_name="GC")
    N50 = models.FloatField(null=True, verbose_name="N50")

    old = models.CharField(max_length=200,null=True, verbose_name="old")


class Gene(models.Model):
    objects = models.Manager()
    MAG = models.ForeignKey(MAG, on_delete=models.CASCADE, verbose_name="MAG", related_name='genes', db_index=True)
    id = models.CharField(max_length=100,  unique=True, primary_key=True, verbose_name="id", db_index=True)
    name = models.CharField(max_length=200, unique=True, verbose_name="name")
    length = models.PositiveIntegerField(verbose_name="length")


class Pathway(models.Model):
    objects = models.Manager()
    MAG = models.ForeignKey(MAG, on_delete=models.CASCADE, verbose_name="MAG", related_name='pathways')
    name = models.CharField(max_length=200, verbose_name="name")
    module_ID = models.CharField(max_length=200, verbose_name="module_ID")
    genes_number = models.PositiveIntegerField(null=True, verbose_name="genes_number")
    steps_number = models.PositiveIntegerField(null=True, verbose_name="steps_number")
    genes_percentage = models.FloatField(null=True, verbose_name="genes_percentage")
    steps_percentage = models.FloatField(null=True, verbose_name="steps_percentage")
    genes = models.CharField(max_length=400,  verbose_name="genes",null=True)
    total_genes = models.PositiveIntegerField(null=True, verbose_name="total_genes")
    total_steps = models.PositiveIntegerField(null=True, verbose_name="total_steps")
""""""


