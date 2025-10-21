import urllib.request
import json
import re
import time

gene_symbol = "MAPT"

print(f"{'='*60}")
print(f"MAPT Geni - ClinVar Uyumlu Analiz")
print(f"{'='*60}\n")

print("ClinVar ile uyumlu MAPT izoformunu çekiliyor...")
print("NOT: ClinVar, 441 amino asitlik beyin izoformunu (P10636-8) kullanır\n")

try:
    protein_id = "P10636-8"
    url = f"https://rest.uniprot.org/uniprotkb/{protein_id}.fasta"
    
    with urllib.request.urlopen(url) as response:
        data = response.read().decode('utf-8')
    
    lines = data.strip().split('\n')
    header = lines[0]
    sequence = ''.join(lines[1:])
    
    if 'sapiens' not in header.lower() and 'human' not in header.lower():
        raise Exception("Yanlış organizma - insan proteini değil!")
    
    print(f"✓ UniProt ID: {protein_id}")
    print(f"✓ İzoform: Tau-F (441 aa) - Beyinde en yaygın izoform")
    print(f"✓ Organizmanın doğrulandı: İnsan (Homo sapiens)")
    print(f"✓ ClinVar uyumlu: Bu izoform ClinVar tarafından kullanılır")

except Exception as e:
    print(f"Hata: {e}")
    print("Alternatif: P10636-8 çekilemedi, temel izoforma geçiliyor...")
    
    protein_id = "P10636"
    url = f"https://rest.uniprot.org/uniprotkb/{protein_id}.fasta"
    
    with urllib.request.urlopen(url) as response:
        data = response.read().decode('utf-8')
    
    lines = data.strip().split('\n')
    header = lines[0]
    sequence = ''.join(lines[1:])
    
    print(f"⚠ UYARI: Kanonik izoform kullanılıyor, pozisyon uyumsuzlukları olabilir")

print(f"\n{'='*60}")
print(f"Protein Bilgileri")
print(f"{'='*60}")
print(f"Protein ID: {protein_id}")
print(f"Başlık: {header}")
print(f"Dizi Uzunluğu: {len(sequence)} amino asit")

if 'sapiens' in header.lower() or 'human' in header.lower():
    print(f"✓ Organizma: Homo sapiens (İnsan) - DOĞRU")
else:
    print(f"⚠ UYARI: Organizma doğrulanamadı!")
    
if len(sequence) == 441:
    print(f"✓ Protein uzunluğu: {len(sequence)} aa - ClinVar referansı ile TAM UYUMLU")
elif len(sequence) >= 352 and len(sequence) <= 758:
    print(f"⚠ Protein uzunluğu: {len(sequence)} aa - MAPT izoformu ama ClinVar uyumsuzluğu olabilir")
    print(f"   ClinVar, 441 aa'lık Tau-F izoformunu kullanır")
else:
    print(f"⚠ UYARI: Beklenmedik protein uzunluğu: {len(sequence)} aa")
    print(f"   MAPT izoformları normalde 352-758 amino asit arasındadır.")
print(f"\nProtein Dizisi:")
print(sequence)

print(f"\nAmino Asit Dağılımı:")
amino_acids = {}
for aa in sequence:
    amino_acids[aa] = amino_acids.get(aa, 0) + 1

for aa in sorted(amino_acids.keys()):
    print(f"{aa}: {amino_acids[aa]} ({amino_acids[aa]/len(sequence)*100:.2f}%)")

print(f"\n{'='*60}")
print(f"ClinVar Varyant Bilgisi ({gene_symbol} geni)")
print(f"{'='*60}")

try:
    esearch_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term={gene_symbol}[gene]&retmode=json&retmax=10000"
    
    with urllib.request.urlopen(esearch_url) as response:
        search_data = json.loads(response.read().decode('utf-8'))
    
    total_variants = int(search_data['esearchresult']['count'])
    print(f"\nToplam Varyant Sayısı: {total_variants}")
    
    pathogenic_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term={gene_symbol}[gene]+AND+(Pathogenic[CLNSIG]+OR+Likely+pathogenic[CLNSIG])&retmode=json&retmax=1"
    
    with urllib.request.urlopen(pathogenic_url) as response:
        pathogenic_data = json.loads(response.read().decode('utf-8'))
    
    pathogenic_count = int(pathogenic_data['esearchresult']['count'])
    
    all_variant_ids = []
    retmax = 500
    for start in range(0, pathogenic_count, retmax):
        fetch_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term={gene_symbol}[gene]+AND+(Pathogenic[CLNSIG]+OR+Likely+pathogenic[CLNSIG])&retmode=json&retstart={start}&retmax={retmax}"
        
        with urllib.request.urlopen(fetch_url) as response:
            fetch_data = json.loads(response.read().decode('utf-8'))
        
        all_variant_ids.extend(fetch_data['esearchresult']['idlist'])
        
        if start + retmax < pathogenic_count:
            time.sleep(0.34)
    
    variant_ids = all_variant_ids
    
    print(f"Patojenik Varyant Sayısı: {pathogenic_count}")
    
    if total_variants > 0:
        print(f"Patojenik Oran: {pathogenic_count/total_variants*100:.2f}%")
    
    print(f"\n{'='*60}")
    print(f"Patojenik Varyantlardaki Amino Asit Değişiklikleri")
    print(f"{'='*60}\n")
    
    amino_acid_changes = {}
    position_mismatches = []
    nonsense_variants = {}
    frameshift_variants = {}
    deletion_variants = {}
    insertion_variants = {}
    duplication_variants = {}
    other_variants = []
    unparsed_variants = []
    non_protein_variants = []
    
    print(f"Toplam {len(variant_ids)} patojenik varyant analiz ediliyor...\n")
    print(f"Referans Protein Uzunluğu: {len(sequence)} amino asit\n")
    
    batch_size = 100
    for i in range(0, len(variant_ids), batch_size):
        batch_ids = variant_ids[i:i+batch_size]
        ids_str = ','.join(batch_ids)
        
        print(f"İşleniyor: {i+1}-{min(i+batch_size, len(variant_ids))} / {len(variant_ids)}", end='\r')
        
        esummary_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=clinvar&id={ids_str}&retmode=json"
        
        with urllib.request.urlopen(esummary_url) as response:
            summary_data = json.loads(response.read().decode('utf-8'))
        
        time.sleep(0.34)
        
        for variant_id in batch_ids:
            if variant_id in summary_data['result']:
                variant_info = summary_data['result'][variant_id]
                title = variant_info.get('title', '')
                
                aa_map = {
                    'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D',
                    'Cys': 'C', 'Gln': 'Q', 'Glu': 'E', 'Gly': 'G',
                    'His': 'H', 'Ile': 'I', 'Leu': 'L', 'Lys': 'K',
                    'Met': 'M', 'Phe': 'F', 'Pro': 'P', 'Ser': 'S',
                    'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V',
                    'Ter': '*', 'X': 'X'
                }
                
                parsed = False
                
                if 'p.' in title:
                    missense_long = re.search(r'p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})', title)
                    missense_short = re.search(r'p\.([A-Z])(\d+)([A-Z])(?![a-z])', title)
                    nonsense_long = re.search(r'p\.([A-Z][a-z]{2})(\d+)(?:Ter|\*)', title)
                    nonsense_short = re.search(r'p\.([A-Z])(\d+)\*', title)
                    frameshift = re.search(r'p\.([A-Z][a-z]{2}|[A-Z])(\d+)(?:[A-Z][a-z]{2}|[A-Z])?fs', title)
                    deletion_single = re.search(r'p\.([A-Z][a-z]{2}|[A-Z])(\d+)del', title)
                    deletion_range = re.search(r'p\.([A-Z][a-z]{2}|[A-Z])(\d+)_([A-Z][a-z]{2}|[A-Z])(\d+)del', title)
                    insertion = re.search(r'p\.([A-Z][a-z]{2}|[A-Z])(\d+)_([A-Z][a-z]{2}|[A-Z])(\d+)ins', title)
                    duplication = re.search(r'p\.([A-Z][a-z]{2}|[A-Z])(\d+)dup', title)
                else:
                    missense_long = missense_short = nonsense_long = nonsense_short = None
                    frameshift = deletion_single = deletion_range = insertion = duplication = None
                
                if missense_long or missense_short:
                    match = missense_long if missense_long else missense_short
                    original = match.group(1)
                    position = match.group(2)
                    mutated = match.group(3)
                    
                    orig_aa = aa_map.get(original, original)
                    mut_aa = aa_map.get(mutated, mutated)
                    pos_int = int(position)
                    
                    if pos_int > len(sequence):
                        position_mismatches.append({
                            'position': pos_int,
                            'change': f"{orig_aa}{position}{mut_aa}",
                            'title': title,
                            'type': 'missense'
                        })
                    elif pos_int <= len(sequence):
                        actual_aa = sequence[pos_int - 1]
                        if actual_aa != orig_aa:
                            position_mismatches.append({
                                'position': pos_int,
                                'change': f"{orig_aa}{position}{mut_aa}",
                                'title': title,
                                'expected': orig_aa,
                                'found': actual_aa,
                                'type': 'missense'
                            })
                        else:
                            key = f"{orig_aa}{position}{mut_aa}"
                            if key not in amino_acid_changes:
                                amino_acid_changes[key] = {
                                    'position': pos_int,
                                    'original': orig_aa,
                                    'mutated': mut_aa,
                                    'count': 0,
                                    'type': 'missense'
                                }
                            amino_acid_changes[key]['count'] += 1
                    parsed = True
                
                elif nonsense_long or nonsense_short:
                    match = nonsense_long if nonsense_long else nonsense_short
                    original = match.group(1)
                    position = match.group(2)
                    
                    orig_aa = aa_map.get(original, original)
                    pos_int = int(position)
                    
                    key = f"{orig_aa}{position}*"
                    if key not in nonsense_variants:
                        nonsense_variants[key] = {
                            'position': pos_int,
                            'original': orig_aa,
                            'count': 0
                        }
                    nonsense_variants[key]['count'] += 1
                    parsed = True
                
                elif frameshift:
                    original = frameshift.group(1)
                    position = frameshift.group(2)
                    
                    orig_aa = aa_map.get(original, original)
                    pos_int = int(position)
                    
                    key = f"{orig_aa}{position}fs"
                    if key not in frameshift_variants:
                        frameshift_variants[key] = {
                            'position': pos_int,
                            'original': orig_aa,
                            'count': 0
                        }
                    frameshift_variants[key]['count'] += 1
                    parsed = True
                
                elif deletion_single:
                    original = deletion_single.group(1)
                    position = deletion_single.group(2)
                    
                    orig_aa = aa_map.get(original, original)
                    pos_int = int(position)
                    
                    key = f"{orig_aa}{position}del"
                    if key not in deletion_variants:
                        deletion_variants[key] = {
                            'position': pos_int,
                            'original': orig_aa,
                            'count': 0,
                            'range': None
                        }
                    deletion_variants[key]['count'] += 1
                    parsed = True
                
                elif deletion_range:
                    start_aa = deletion_range.group(1)
                    start_pos = deletion_range.group(2)
                    end_aa = deletion_range.group(3)
                    end_pos = deletion_range.group(4)
                    
                    start_aa_short = aa_map.get(start_aa, start_aa)
                    end_aa_short = aa_map.get(end_aa, end_aa)
                    
                    key = f"{start_aa_short}{start_pos}_{end_aa_short}{end_pos}del"
                    if key not in deletion_variants:
                        deletion_variants[key] = {
                            'position': int(start_pos),
                            'original': f"{start_aa_short}-{end_aa_short}",
                            'count': 0,
                            'range': (int(start_pos), int(end_pos))
                        }
                    deletion_variants[key]['count'] += 1
                    parsed = True
                
                elif insertion:
                    start_aa = insertion.group(1)
                    start_pos = insertion.group(2)
                    end_aa = insertion.group(3)
                    end_pos = insertion.group(4)
                    
                    start_aa_short = aa_map.get(start_aa, start_aa)
                    end_aa_short = aa_map.get(end_aa, end_aa)
                    
                    key = f"{start_aa_short}{start_pos}_{end_aa_short}{end_pos}ins"
                    if key not in insertion_variants:
                        insertion_variants[key] = {
                            'position': int(start_pos),
                            'count': 0
                        }
                    insertion_variants[key]['count'] += 1
                    parsed = True
                
                elif duplication:
                    original = duplication.group(1)
                    position = duplication.group(2)
                    
                    orig_aa = aa_map.get(original, original)
                    pos_int = int(position)
                    
                    key = f"{orig_aa}{position}dup"
                    if key not in duplication_variants:
                        duplication_variants[key] = {
                            'position': pos_int,
                            'original': orig_aa,
                            'count': 0
                        }
                    duplication_variants[key]['count'] += 1
                    parsed = True
                
                if not parsed and 'p.' in title:
                    unparsed_variants.append(title)
                elif not parsed and 'p.' not in title:
                    non_protein_variants.append(title)
    
    sorted_changes = sorted(amino_acid_changes.items(), key=lambda x: x[1]['position'])
    sorted_nonsense = sorted(nonsense_variants.items(), key=lambda x: x[1]['position'])
    sorted_frameshift = sorted(frameshift_variants.items(), key=lambda x: x[1]['position'])
    sorted_deletion = sorted(deletion_variants.items(), key=lambda x: x[1]['position'])
    sorted_insertion = sorted(insertion_variants.items(), key=lambda x: x[1]['position'])
    sorted_duplication = sorted(duplication_variants.items(), key=lambda x: x[1]['position'])
    
    total_parsed = (len(amino_acid_changes) + len(nonsense_variants) + 
                    len(frameshift_variants) + len(deletion_variants) + 
                    len(insertion_variants) + len(duplication_variants))
    
    total_accounted = total_parsed + len(position_mismatches) + len(unparsed_variants) + len(non_protein_variants)
    
    print(f"\n\n{'='*60}")
    print(f"Analiz Sonuçları - Mutasyon Türlerine Göre")
    print(f"{'='*60}")
    print(f"Toplam Patojenik Varyant: {pathogenic_count}")
    print(f"Çekilen Varyant ID: {len(variant_ids)}")
    print(f"\nMutasyon Türleri:")
    print(f"  - Missense (Yanlış Anlamlı): {len(amino_acid_changes)}")
    print(f"  - Nonsense (Anlamsız/Stop): {len(nonsense_variants)}")
    print(f"  - Frameshift (Çerçeve Kayması): {len(frameshift_variants)}")
    print(f"  - Deletion (Delesyon): {len(deletion_variants)}")
    print(f"  - Insertion (İnsersiyon): {len(insertion_variants)}")
    print(f"  - Duplication (Duplikasyon): {len(duplication_variants)}")
    print(f"\nToplam Parse Edilen: {total_parsed}")
    print(f"İzoform Uyumsuzluğu: {len(position_mismatches)}")
    print(f"Parse Edilemeyen (p. notasyonu ile): {len(unparsed_variants)}")
    print(f"Protein Dışı Varyantlar (p. notasyonu yok): {len(non_protein_variants)}")
    print(f"\nToplam Hesaba Katılan: {total_accounted}")
    print(f"Eksik Varyant: {pathogenic_count - total_accounted}")
    print(f"{'='*60}\n")
    
    if position_mismatches:
        print(f"\n{'='*60}")
        print(f"Bilgi: Farklı İzoformlara Ait Varyantlar Tespit Edildi")
        print(f"{'='*60}")
        print(f"\nAşağıdaki {len(position_mismatches)} varyant, mevcut protein izoformu")
        print(f"(UniProt: {protein_id}, {len(sequence)} aa) ile uyumlu değil.\n")
        print(f"Bu varyantlar muhtemelen MAPT'nin diğer izoformlarına aittir:")
        
        out_of_range = [m for m in position_mismatches if m['position'] > len(sequence)]
        aa_mismatch = [m for m in position_mismatches if m['position'] <= len(sequence)]
        
        if out_of_range:
            print(f"\n1) Daha Uzun İzoformlara Ait ({len(out_of_range)} varyant):")
            print(f"   (Pozisyon > {len(sequence)}, muhtemelen 758 aa izoformda)")
            for m in out_of_range[:10]:
                print(f"   - {m['change']} (Pozisyon {m['position']})")
            if len(out_of_range) > 10:
                print(f"   ... ve {len(out_of_range) - 10} varyant daha")
        
        if aa_mismatch:
            print(f"\n2) Alternatif Splicing Bölgelerindeki Varyantlar ({len(aa_mismatch)} varyant):")
            print(f"   (Amino asit uyuşmazlığı - farklı ekson kullanımı)")
            for m in aa_mismatch[:10]:
                if 'expected' in m and 'found' in m:
                    print(f"   - {m['change']}: Pos {m['position']} - Beklenen '{m['expected']}', Bu izoformda '{m['found']}'")
                else:
                    print(f"   - {m['change']}: Pos {m['position']}")
            if len(aa_mismatch) > 10:
                print(f"   ... ve {len(aa_mismatch) - 10} varyant daha")
        
        print(f"\n{'='*60}")
        print(f"NOT: MAPT geninin 6 farklı izoformu vardır (352-758 aa).")
        print(f"ClinVar varyantları farklı izoformlar için raporlanmış olabilir.")
        print(f"Şu anki analiz {protein_id} ({len(sequence)} aa) izoformu içindir.")
        print(f"{'='*60}\n")
    
    if sorted_changes:
        print(f"\n{'='*60}")
        print(f"1. Missense Mutasyonlar ({len(sorted_changes)})")
        print(f"{'='*60}\n")
        for change, info in sorted_changes:
            print(f"Pozisyon {info['position']:4d}: {info['original']} → {info['mutated']}  ({change})")
        
        print(f"\nEn Sık Görülen Missense Mutasyonlar:")
        top_missense = sorted(amino_acid_changes.items(), key=lambda x: x[1]['count'], reverse=True)[:5]
        for change, info in top_missense:
            print(f"  {change}: {info['count']} kez")
    
    if sorted_nonsense:
        print(f"\n{'='*60}")
        print(f"2. Nonsense Mutasyonlar (Stop Codon) ({len(sorted_nonsense)})")
        print(f"{'='*60}\n")
        for change, info in sorted_nonsense:
            print(f"Pozisyon {info['position']:4d}: {info['original']} → STOP  ({change})")
    
    if sorted_frameshift:
        print(f"\n{'='*60}")
        print(f"3. Frameshift Mutasyonlar ({len(sorted_frameshift)})")
        print(f"{'='*60}\n")
        for change, info in sorted_frameshift:
            print(f"Pozisyon {info['position']:4d}: {info['original']} frameshift  ({change})")
    
    if sorted_deletion:
        print(f"\n{'='*60}")
        print(f"4. Delesyon Mutasyonları ({len(sorted_deletion)})")
        print(f"{'='*60}\n")
        for change, info in sorted_deletion:
            if info['range']:
                print(f"Pozisyon {info['position']:4d}-{info['range'][1]}: Delesyon  ({change})")
            else:
                print(f"Pozisyon {info['position']:4d}: {info['original']} delesyon  ({change})")
    
    if sorted_insertion:
        print(f"\n{'='*60}")
        print(f"5. İnsersiyon Mutasyonları ({len(sorted_insertion)})")
        print(f"{'='*60}\n")
        for change, info in sorted_insertion:
            print(f"Pozisyon {info['position']:4d}: İnsersiyon  ({change})")
    
    if sorted_duplication:
        print(f"\n{'='*60}")
        print(f"6. Duplikasyon Mutasyonları ({len(sorted_duplication)})")
        print(f"{'='*60}\n")
        for change, info in sorted_duplication:
            print(f"Pozisyon {info['position']:4d}: {info['original']} duplikasyon  ({change})")
    
    if unparsed_variants:
        print(f"\n{'='*60}")
        print(f"UYARI: Parse Edilemeyen Protein Varyantları ({len(unparsed_variants)})")
        print(f"{'='*60}\n")
        print("Bu varyantlar 'p.' notasyonu içeriyor ancak mevcut regex'lerle parse edilemedi:")
        for title in unparsed_variants[:20]:
            print(f"  - {title}")
        if len(unparsed_variants) > 20:
            print(f"  ... ve {len(unparsed_variants) - 20} varyant daha")
    
    if non_protein_variants:
        print(f"\n{'='*60}")
        print(f"Protein Dışı Varyantlar ({len(non_protein_variants)})")
        print(f"{'='*60}\n")
        print("Bu varyantlar 'p.' notasyonu içermiyor (muhtemelen DNA seviyesinde varyantlar):")
        for title in non_protein_variants[:20]:
            print(f"  - {title}")
        if len(non_protein_variants) > 20:
            print(f"  ... ve {len(non_protein_variants) - 20} varyant daha")
    
    print(f"\n{'='*60}")
    print(f"Mutasyon Şiddeti Değerlendirmesi")
    print(f"{'='*60}\n")
    
    total_loss_of_function = len(nonsense_variants) + len(frameshift_variants) + len(deletion_variants)
    total_variants_analyzed = total_parsed
    
    print(f"Protein Fonksiyon Kaybına Neden Olan Mutasyonlar:")
    print(f"  - Nonsense + Frameshift + Delesyon: {total_loss_of_function}")
    print(f"  - Toplam Mutasyonların %{total_loss_of_function/total_variants_analyzed*100:.1f}'i")
    print(f"\nProtein Yapısını Değiştiren Mutasyonlar:")
    print(f"  - Missense: {len(amino_acid_changes)}")
    print(f"  - Toplam Mutasyonların %{len(amino_acid_changes)/total_variants_analyzed*100:.1f}'i")
    print(f"\nNOT: Nonsense, Frameshift ve Delesyon mutasyonları genellikle")
    print(f"Missense mutasyonlardan daha şiddetli patolojik etkilere sahiptir.")
    print(f"{'='*60}")

except Exception as e:
    print(f"ClinVar sorgusu sırasında hata: {str(e)}")

