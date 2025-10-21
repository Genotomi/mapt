# MAPT İzoform-Özgül Varyant Analizcisi

Bu bir Python betiğidir.

## Kodun Amacı

Bu betik, MAPT (Tau) genine ait patojenik varyantları ClinVar veritabanından almak ve bu varyantların **spesifik bir protein izoformu** (beynin en yaygın 441 amino asitlik Tau-F, P10636-8) ile uyumluluğunu doğrulamak için tasarlanmıştır.

Temel amacı, bir varyantın bu spesifik izoformla uyumlu olup olmadığını belirlemek için iki aşamalı bir filtreleme yapmaktır:

1.  Varyantın pozisyonunun 441 amino asitlik dizi dışında olup olmadığını kontrol eder.
2.  Pozisyon doğru olsa bile, amino asidin (örn. S226Y) referans izoformdaki (V226) amino asitle eşleşip eşleşmediğini kontrol eder.

Bu filtreleme, "protein-dışı" varyantları (`p.` notasyonu olmayanlar) ve izoformla uyumsuz olanları ayırarak, yalnızca P10636-8 izoformuyla tam uyumlu olan varyantlardan oluşan temiz bir veri seti sağlar.
