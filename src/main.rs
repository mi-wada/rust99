use std::fmt;

fn main() {}

#[cfg(test)]
mod tests {
    use super::*;
    use std::{collections::HashSet, hash::Hash, iter::FromIterator};

    use rand::{prelude::SliceRandom, thread_rng};

    #[test]
    fn p01() {
        let li = vec![1, 1, 2, 3, 5, 8];
        fn last(v: &Vec<i32>) -> Option<i32> {
            v.last().copied()
        }
        assert_eq!(Some(8), last(&li));
    }

    #[test]
    fn p02() {
        let li = vec![1, 1, 2, 3, 5, 8];
        fn penultimate(v: &Vec<i32>) -> Option<i32> {
            v[..v.len() - 1].last().copied()
        }
        assert_eq!(Some(5), penultimate(&li));
    }

    #[test]
    fn p03() {
        let li = vec![1, 1, 2, 3, 5, 8];
        fn nth(v: &Vec<i32>, n: usize) -> Option<i32> {
            v.iter().nth(n).copied()
        }
        assert_eq!(Some(2), nth(&li, 2));
    }

    #[test]
    fn p04() {
        let li = vec![1, 1, 2, 3, 5, 8];
        fn length(v: &Vec<i32>) -> usize {
            v.len()
        }
        assert_eq!(6, length(&li));
    }

    #[test]
    fn p05() {
        let li = vec![1, 1, 2, 3, 5, 8];
        fn reverse(v: &Vec<i32>) -> Vec<i32> {
            let mut res = Vec::with_capacity(v.len());
            let mut copied_v = v.clone();
            while let Some(x) = copied_v.pop() {
                res.push(x);
            }
            res
        }
        assert_eq!(vec![8, 5, 3, 2, 1, 1], reverse(&li));
    }

    #[test]
    fn p06() {
        fn is_palindrome(v: &Vec<i32>) -> bool {
            for i in 0..v.len() / 2 {
                if v[i] != v[v.len() - 1 - i] {
                    return false;
                }
            }
            true
        }
        let palindrome = vec![1, 2, 3, 2, 1];
        assert_eq!(true, is_palindrome(&palindrome));
        let not_palindrome = vec![1, 3, 3, 2, 1];
        assert_eq!(false, is_palindrome(&not_palindrome));
    }

    #[test]
    fn p07() {
        #[derive(PartialEq, Debug)]
        struct List<T: Copy>(Vec<ListItem<T>>);

        impl<T: Copy> List<T> {
            fn new() -> List<T> {
                List(vec![])
            }

            fn from(items: Vec<ListItem<T>>) -> List<T> {
                List(items)
            }

            fn push(&mut self, element: ListItem<T>) {
                self.0.push(element)
            }

            fn append(&mut self, other: List<T>) {
                for x in other.0.into_iter() {
                    self.0.push(x);
                }
            }
        }

        #[derive(Debug, PartialEq)]
        enum ListItem<T: Copy> {
            Item(T),
            NestedItem(Box<List<T>>),
        }

        fn item<T: Copy>(x: T) -> ListItem<T> {
            ListItem::Item(x)
        }

        fn nested_item<T: Copy>(li: List<T>) -> ListItem<T> {
            ListItem::NestedItem(Box::new(li))
        }

        let not_flat_list = List::from(vec![
            item(10),
            nested_item(List::from(vec![item(20), item(30)])),
        ]);

        fn flatten<T: Copy>(li: List<T>) -> List<T> {
            let mut res = List::new();
            for x in li.0 {
                match x {
                    ListItem::Item(x) => res.push(item(x)),
                    ListItem::NestedItem(li) => res.append(flatten(*li)),
                }
            }
            res
        }

        assert_eq!(
            flatten(not_flat_list),
            List::from(vec![item(10), item(20), item(30)])
        );
    }

    #[test]
    fn p08() {
        let li = vec![
            'a', 'a', 'a', 'a', 'b', 'c', 'c', 'a', 'a', 'd', 'e', 'e', 'e', 'e',
        ];
        fn compress<T>(v: &Vec<T>) -> Vec<T>
        where
            T: Copy + Eq + Hash,
        {
            let mut repeated_detecter = HashSet::<T>::new();
            let mut compressd = vec![];
            for x in v {
                if !repeated_detecter.contains(x) {
                    compressd.push(*x);
                }
                repeated_detecter.clear();
                repeated_detecter.insert(*x);
            }
            compressd
        }
        assert_eq!(compress(&li), vec!['a', 'b', 'c', 'a', 'd', 'e']);
    }

    #[test]
    fn p09() {
        let li = vec![
            'a', 'a', 'a', 'a', 'b', 'c', 'c', 'a', 'a', 'd', 'e', 'e', 'e', 'e',
        ];
        fn pack<T: Copy + PartialEq>(v: &Vec<T>) -> Vec<Vec<T>> {
            let mut packed = vec![vec![]];
            if v.is_empty() {
                return packed;
            }
            let mut prev_value = &v[0];
            for x in v {
                if x == prev_value {
                    packed.last_mut().unwrap().push(*x);
                } else {
                    packed.push(vec![*x]);
                }
                prev_value = x;
            }
            packed
        }
        assert_eq!(
            pack(&li),
            vec![
                vec!['a', 'a', 'a', 'a'],
                vec!['b'],
                vec!['c', 'c'],
                vec!['a', 'a'],
                vec!['d'],
                vec!['e', 'e', 'e', 'e']
            ]
        );
    }

    #[test]
    fn p10() {
        let li = vec![
            'a', 'a', 'a', 'a', 'b', 'c', 'c', 'a', 'a', 'd', 'e', 'e', 'e', 'e',
        ];
        fn pack<T: Copy + PartialEq>(v: &Vec<T>) -> Vec<Vec<T>> {
            let mut packed = vec![vec![]];
            if v.is_empty() {
                return packed;
            }
            let mut prev_value = &v[0];
            for x in v {
                if x == prev_value {
                    packed.last_mut().unwrap().push(*x);
                } else {
                    packed.push(vec![*x]);
                }
                prev_value = x;
            }
            packed
        }
        fn encode<T: Copy + PartialEq>(v: &Vec<T>) -> Vec<(usize, T)> {
            let mut encoded = vec![];
            for x in pack(v) {
                encoded.push((x.len(), *x.last().unwrap()));
            }
            encoded
        }
        assert_eq!(
            encode(&li),
            vec![(4, 'a'), (1, 'b'), (2, 'c'), (2, 'a'), (1, 'd'), (4, 'e')]
        );
    }

    fn p12() {
        let li = vec![(4, 'a'), (1, 'b'), (2, 'c'), (2, 'a'), (1, 'd'), (4, 'e')];
        fn decode<T: Copy + PartialEq>(v: &Vec<(usize, T)>) -> Vec<T> {
            let mut decoded = vec![];
            for x in v {
                for _ in 0..x.0 {
                    decoded.push(x.1);
                }
            }
            decoded
        }
        assert_eq!(
            decode(&li),
            vec!['a', 'a', 'a', 'a', 'b', 'c', 'c', 'a', 'a', 'd', 'e', 'e', 'e', 'e']
        );
    }

    #[test]
    fn p13() {
        let li = vec![
            'a', 'a', 'a', 'a', 'b', 'c', 'c', 'a', 'a', 'd', 'e', 'e', 'e', 'e',
        ];
        fn encode<T: Copy + PartialEq>(v: &Vec<T>) -> Vec<(usize, T)> {
            if v.is_empty() {
                vec![]
            } else {
                let mut encoded = vec![(1, v[0].clone())];
                for x in v[1..].iter() {
                    if x == &encoded.last().unwrap().1 {
                        encoded.last_mut().unwrap().0 += 1;
                    } else {
                        encoded.push((1, *x));
                    }
                }
                encoded
            }
        }
        assert_eq!(
            encode(&li),
            [(4, 'a'), (1, 'b'), (2, 'c'), (2, 'a'), (1, 'd'), (4, 'e')]
        );
    }

    #[test]
    fn p14() {
        let li = vec!['a', 'b', 'c', 'c', 'd'];
        fn duplicate<T>(li: &Vec<T>) -> Vec<&T> {
            let mut duplicated = vec![];
            let duplicated_count = 2;
            for x in li {
                for _ in 0..duplicated_count {
                    duplicated.push(x);
                }
            }
            duplicated
        }
        assert_eq!(
            duplicate(&li),
            vec![&'a', &'a', &'b', &'b', &'c', &'c', &'c', &'c', &'d', &'d']
        );
    }

    #[test]
    fn p15() {
        let li = vec!['a', 'b', 'c', 'c', 'd'];
        fn duplicate<T>(li: &Vec<T>, duplicated_count: usize) -> Vec<&T> {
            let mut duplicated = vec![];
            for x in li {
                for _ in 0..duplicated_count {
                    duplicated.push(x);
                }
            }
            duplicated
        }
        assert_eq!(
            duplicate(&li, 2),
            vec![&'a', &'a', &'b', &'b', &'c', &'c', &'c', &'c', &'d', &'d']
        );
    }

    #[test]
    fn p16() {
        let li = vec!['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k'];
        fn drop_every_nth<T>(li: &Vec<T>, n: usize) -> Vec<&T> {
            let mut droped = vec![];
            for (i, x) in li.iter().enumerate() {
                if (i + 1) % n != 0 {
                    droped.push(x);
                }
            }
            droped
        }
        assert_eq!(
            drop_every_nth(&li, 3),
            vec![&'a', &'b', &'d', &'e', &'g', &'h', &'j', &'k']
        );
    }

    #[test]
    fn p17() {
        let li = vec!['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k'];
        fn split<T>(v: &Vec<T>, splited_at: usize) -> (&[T], &[T]) {
            v.split_at(splited_at)
        }
        assert_eq!(
            split(&li, 3),
            (
                vec!['a', 'b', 'c'].as_slice(),
                vec!['d', 'e', 'f', 'g', 'h', 'i', 'j', 'k'].as_slice()
            )
        )
    }

    #[test]
    fn p18() {
        let li = vec!['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k'];
        fn slice<T>(v: &Vec<T>, begin: usize, end: usize) -> &[T] {
            &v[begin..end]
        }
        assert_eq!(slice(&li, 3, 7), ['d', 'e', 'f', 'g']);
    }

    #[test]
    // DIFFICULT
    fn p19() {
        let li = vec!['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k'];
        fn rotate<T: Clone>(v: &Vec<T>, from: isize) -> Vec<T> {
            let from = if from < 0 {
                v.len() - (from - 1) as usize % v.len()
            } else {
                from as usize % v.len()
            };
            let mut rotated = v[from..v.len()].to_vec();
            rotated.append(&mut v[0..from].to_vec());
            rotated
        }
        assert_eq!(
            rotate(&li, -2),
            vec!['j', 'k', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i']
        );
        assert_eq!(
            rotate(&li, 3),
            vec!['d', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'a', 'b', 'c']
        );
    }

    #[test]
    fn p20() {
        let mut li = vec!['a', 'b', 'c', 'd'];
        fn remove_at<T: Copy>(v: &mut Vec<T>, at: usize) -> (Vec<T>, T) {
            let removed_item = v.remove(at);
            (v.clone(), removed_item)
        }
        assert_eq!(remove_at(&mut li, 1), (vec!['a', 'c', 'd'], 'b'));
    }

    #[test]
    fn p21() {
        let mut li = vec!['a', 'b', 'c', 'd'];
        fn insert_at<T: Clone>(v: &mut Vec<T>, at: usize, item: T) -> Vec<T> {
            v.insert(at, item);
            v.clone()
        }
        assert_eq!(insert_at(&mut li, 1, 'n'), vec!['a', 'n', 'b', 'c', 'd']);
    }

    #[test]
    fn p22() {
        fn range(begin: isize, end: isize) -> Vec<isize> {
            (begin..=end).collect()
        }
        assert_eq!(range(0, 3), vec![0, 1, 2, 3]);
    }

    #[test]
    fn p23() {
        let li = vec!['j', 'k', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i'];
        fn random_select<T>(v: &Vec<T>, len: usize) -> Vec<&T> {
            let mut random_selected = vec![];
            let mut rng = thread_rng();
            for _ in 0..len {
                random_selected.push(v.choose(&mut rng).unwrap());
            }
            random_selected
        }
        let random_selected = random_select(&li, 3);
        assert_eq!(random_selected.len(), 3);
        assert_eq!(
            HashSet::<&char>::from_iter(random_selected.into_iter())
                .is_subset(&HashSet::<&char>::from_iter(li.iter())),
            true
        );
    }

    #[test]
    fn p24() {
        fn lotto(len: usize, end: usize) -> Vec<usize> {
            let mut res = vec![];
            let candidates = (1..end).collect::<Vec<usize>>();
            let mut rng = thread_rng();
            for _ in 0..len {
                res.push(*candidates.choose(&mut rng).unwrap());
            }
            res
        }

        let (len, end) = (6, 49);
        assert_eq!(lotto(len, end).len(), len);
        assert_eq!(
            HashSet::<usize>::from_iter(lotto(len, end).into_iter())
                .is_subset(&HashSet::from_iter(1..end)),
            true
        );
    }

    #[test]
    fn p25() {
        let li = vec!['a', 'b', 'c', 'd', 'e', 'f'];
        fn random_permute<T: Clone>(v: &Vec<T>) -> Vec<T> {
            let mut rng = thread_rng();
            let mut cloned_v = v.clone();
            cloned_v.shuffle(&mut rng);
            cloned_v
        }
        assert_eq!(
            HashSet::<&char>::from_iter(random_permute(&li).iter())
                == HashSet::from_iter(li.iter()),
            true
        );
    }

    #[test]
    fn p26() {
        let li = vec!['a', 'b', 'c', 'd', 'e'];
        fn combinations<T: Copy>(k: usize, li: &[T]) -> Vec<Vec<T>> {
            fn comb_rec<T: Copy>(k: usize, remaining: &[T], acc: &Vec<T>) -> Vec<Vec<T>> {
                if k == 1 {
                    remaining
                        .into_iter()
                        .map(|&x| {
                            let mut acc2 = acc.clone();
                            acc2.push(x);
                            acc2
                        })
                        .collect()
                } else {
                    let mut res = vec![];
                    for i in 0..remaining.len() - k + 1 {
                        let (left, right) = remaining.split_at(i + 1);
                        let x = left[left.len() - 1];
                        let mut acc2 = acc.clone();
                        acc2.push(x);
                        res.extend(comb_rec(k - 1, &right, &acc2));
                    }
                    res
                }
            }
            comb_rec(k, &li, &vec![])
        }
    }

    #[test]
    fn p27() {
        let li = vec![
            "Aldo", "Beat", "Carla", "David", "Evi", "Flip", "Gary", "Hugo", "Ida",
        ];
        fn combinations<T: Copy>(k: usize, li: &[T]) -> Vec<Vec<T>> {
            fn comb_rec<T: Copy>(k: usize, remaining: &[T], acc: &Vec<T>) -> Vec<Vec<T>> {
                if k == 1 {
                    remaining
                        .into_iter()
                        .map(|&x| {
                            let mut acc2 = acc.clone();
                            acc2.push(x);
                            acc2
                        })
                        .collect()
                } else {
                    let mut res = vec![];
                    for i in 0..remaining.len() - k + 1 {
                        let (left, right) = remaining.split_at(i + 1);
                        let x = left[left.len() - 1];
                        let mut acc2 = acc.clone();
                        acc2.push(x);
                        res.extend(comb_rec(k - 1, &right, &acc2));
                    }
                    res
                }
            }
            comb_rec(k, &li, &vec![])
        }
        fn group3<T: Copy>(v: &Vec<T>) -> Vec<Vec<Vec<T>>> {
            let mut res = vec![];
            let combs = combinations(v.len(), v);
            for comb in combs {
                res.push(vec![
                    comb[0..2].to_vec(),
                    comb[2..5].to_vec(),
                    comb[5..9].to_vec(),
                ]);
            }
            res
        }
        let groups = group3(&li);
    }

    #[test]
    fn p28() {
        let li = vec![
            vec!['a', 'b', 'c'],
            vec!['d', 'e'],
            vec!['f', 'g', 'h'],
            vec!['d', 'e'],
            vec!['i', 'j', 'k', 'l'],
            vec!['m', 'n'],
            vec!['o'],
        ];

        fn lsort<T: Clone>(v: &Vec<Vec<T>>) -> Vec<Vec<T>> {
            let mut res = v.clone();
            res.sort_by(|a, b| a.len().partial_cmp(&b.len()).unwrap());
            res
        }

        assert_eq!(
            lsort(&li),
            vec![
                vec!['o'],
                vec!['d', 'e'],
                vec!['d', 'e'],
                vec!['m', 'n'],
                vec!['a', 'b', 'c'],
                vec!['f', 'g', 'h'],
                vec!['i', 'j', 'k', 'l']
            ]
        );
    }

    fn p29() {
        struct PrimeJudgement {
            max_num: usize,
            primes: HashSet<usize>,
        }

        impl PrimeJudgement {
            fn new() -> Self {
                Self {
                    max_num: 0,
                    primes: HashSet::<usize>::new(),
                }
            }
            fn is_prime(&mut self, n: usize) -> bool {
                if n > self.max_num {
                    for x in self.max_num + 1..n {
                        let mut x_is_prime = true;
                        for prime in self.primes.iter() {
                            if x % *prime == 0 {
                                x_is_prime = false;
                                break;
                            }
                        }
                        if x_is_prime {
                            self.primes.insert(x);
                        }
                    }
                }
                self.primes.contains(&n)
            }
        }

        let mut prime_judgement = PrimeJudgement::new();
        assert_eq!(prime_judgement.is_prime(53), true);
        assert_eq!(prime_judgement.is_prime(1957), false);
        assert_eq!(prime_judgement.is_prime(3571), true);
    }

    #[test]
    fn p32() {
        fn gcd(a: usize, b: usize) -> usize {
            if b == 0 {
                a
            } else {
                gcd(b, a % b)
            }
        }
        assert_eq!(gcd(36, 63), 9);
    }

    #[test]
    fn p33() {
        fn gcd(a: usize, b: usize) -> usize {
            if b == 0 {
                a
            } else {
                gcd(b, a % b)
            }
        }
        fn is_coprime_to(a: usize, b: usize) -> bool {
            gcd(a, b) == 1
        }

        assert_eq!(is_coprime_to(35, 64), true);
        assert_eq!(is_coprime_to(35, 65), false);
    }

    #[test]
    fn p34() {
        fn gcd(a: usize, b: usize) -> usize {
            if b == 0 {
                a
            } else {
                gcd(b, a % b)
            }
        }
        fn is_coprime_to(a: usize, b: usize) -> bool {
            gcd(a, b) == 1
        }
        fn totient(n: usize) -> usize {
            (1..=n).filter(|x| is_coprime_to(n, *x)).count()
        }

        assert_eq!(totient(10), 4);
    }

    // 面倒なのでskip
    // next: p46
    #[test]
    fn p49() {
        fn to_binary_str(n: usize, digit: usize) -> String {
            let mut n = n;
            let mut rev_res = "".to_string();
            while rev_res.len() < digit {
                rev_res.push_str(&(n % 2).to_string());
                n /= 2;
            }
            let res = {
                let mut res = "".to_string();
                while let Some(c) = rev_res.pop() {
                    res.push(c);
                }
                res
            };
            res
        }
        fn gray(n: usize) -> Vec<String> {
            (0..(2i32.pow(n as u32)))
                .map(|x| {
                    let code = (x >> 1) ^ x;
                    to_binary_str(code as usize, n)
                })
                .collect()
        }

        assert_eq!(
            gray(3),
            vec![
                "000".to_string(),
                "001".to_string(),
                "011".to_string(),
                "010".to_string(),
                "110".to_string(),
                "111".to_string(),
                "101".to_string(),
                "100".to_string()
            ]
        )
    }

    #[test]
    fn p50() {
        enum Node {
            Leaf(usize, char),
            Internal(usize, Box<Node>, Box<Node>),
        }

        impl Node {
            fn value(&self) -> usize {
                match *self {
                    Node::Internal(v, _, _) => v,
                    Node::Leaf(v, _) => v,
                }
            }
        }

        fn huffman(origin: &Vec<(char, usize)>) -> Vec<(char, String)> {
            let mut nodes = origin
                .iter()
                .map(|(c, v)| Node::Leaf(*v, *c))
                .collect::<Vec<Node>>();
            while nodes.len() > 1 {
                nodes.sort_by(|a, b| b.value().cmp(&a.value()));
                let (smallest, second_smallest) = (nodes.pop().unwrap(), nodes.pop().unwrap());

                nodes.push(Node::Internal(
                    smallest.value() + second_smallest.value(),
                    Box::new(smallest),
                    Box::new(second_smallest),
                ));
            }
            let huffman_tree = &nodes[0];
            huffman_tbl(huffman_tree, "".to_string())
        }

        fn huffman_tbl(huffman_tree: &Node, acc: String) -> Vec<(char, String)> {
            match huffman_tree {
                Node::Leaf(_, c) => {
                    vec![(*c, acc)]
                }
                Node::Internal(_, l, r) => {
                    let mut res = vec![];
                    res.extend(huffman_tbl(l, {
                        let mut acc_l = acc.clone();
                        acc_l.push('0');
                        acc_l
                    }));
                    res.extend(huffman_tbl(r, {
                        let mut acc_r = acc.clone();
                        acc_r.push('1');
                        acc_r
                    }));
                    res
                }
            }
        }

        let symbols = vec![
            ('a', 45),
            ('b', 13),
            ('c', 12),
            ('d', 16),
            ('e', 9),
            ('f', 5),
        ];
        assert_eq!(
            huffman(&symbols),
            vec![
                ('a', "0".to_string()),
                ('c', "100".to_string()),
                ('b', "101".to_string()),
                ('f', "1100".to_string()),
                ('e', "1101".to_string()),
                ('d', "111".to_string()),
            ]
        )
    }

    #[test]
    fn p55() {
        fn cbal_trees(n: usize, symbol: char) -> Vec<Tree<char>> {
            if n == 0 {
                vec![Tree::end()]
            } else if (n - 1) % 2 == 0 {
                let subtrees = cbal_trees((n - 1) / 2, symbol);

                subtrees
                    .iter()
                    .flat_map(|s| {
                        subtrees
                            .iter()
                            .map(move |s_| Tree::node(symbol, s.clone(), s_.clone()))
                    })
                    .collect()
            } else {
                let subtrees1 = cbal_trees(n / 2, symbol);
                let subtrees2 = cbal_trees(n / 2 - 1, symbol);

                subtrees1
                    .iter()
                    .flat_map(|s1| {
                        subtrees2.iter().flat_map(move |s2| {
                            vec![
                                Tree::node(symbol, s1.clone(), s2.clone()),
                                Tree::node(symbol, s2.clone(), s1.clone()),
                            ]
                        })
                    })
                    .collect()
            }
        }

        // println!("{:#?}", cbal_trees(4, 'x'));
    }

    #[test]
    fn p56() {
        // not interested in value
        fn is_symmetric<T: std::fmt::Display>(tree: &Tree<T>) -> bool {
            match tree {
                Tree::End => true,
                Tree::Node {
                    value: _,
                    left,
                    right,
                } => left.is_mirror_of(right),
            }
        }
        let symmetric_tree = Tree::node('x', Tree::leaf('y'), Tree::leaf('z'));
        assert_eq!(is_symmetric(&symmetric_tree), true);
        let symmetric_tree = Tree::node('x', Tree::leaf('x'), Tree::end());
        assert_eq!(is_symmetric(&symmetric_tree), false);
    }

    #[test]
    fn p57() {
        fn add_value<T: std::fmt::Display + std::cmp::PartialOrd + Copy>(
            tree: &mut Tree<T>,
            value: T,
        ) {
            match tree {
                Tree::End => *tree = Tree::leaf(value),
                Tree::Node {
                    value: tree_value,
                    left,
                    right,
                } => {
                    if value < *tree_value {
                        add_value(left, value);
                    } else {
                        add_value(right, value);
                    }
                }
            }
        }

        fn from_list<T: std::fmt::Display + Copy + std::cmp::PartialOrd>(v: &Vec<T>) -> Tree<T> {
            let mut tree = Tree::end();
            for x in v {
                add_value(&mut tree, *x);
            }
            tree
        }

        let mut tree = Tree::end();
        add_value(&mut tree, 2);
        add_value(&mut tree, 3);
        add_value(&mut tree, 0);
        assert_eq!(tree, Tree::node(2, Tree::leaf(0), Tree::leaf(3)));
    }

    #[test]
    fn p58() {
        fn balanced_trees<T: fmt::Display + Copy>(n: usize, value: T) -> Vec<Tree<T>> {
            if n == 0 {
                vec![Tree::End]
            } else if (n - 1) % 2 == 0 {
                let balanced_subtrees = balanced_trees((n - 1) / 2, value);
                balanced_subtrees
                    .iter()
                    .flat_map(|s| {
                        balanced_subtrees
                            .iter()
                            .map(move |s_| Tree::node(value, s.clone(), s_.clone()))
                    })
                    .collect()
            } else {
                let balanced_subtrees1 = balanced_trees(n / 2, value);
                let balanced_subtrees2 = balanced_trees(n / 2 - 1, value);
                balanced_subtrees1
                    .iter()
                    .flat_map(|s1| {
                        balanced_subtrees2.iter().flat_map(move |s2| {
                            vec![
                                Tree::node(value, s1.clone(), s2.clone()),
                                Tree::node(value, s2.clone(), s1.clone()),
                            ]
                        })
                    })
                    .collect()
            }
        }

        fn is_symmetric<T: fmt::Display>(tree: &Tree<T>) -> bool {
            match tree {
                Tree::End => true,
                Tree::Node {
                    value: _,
                    left,
                    right,
                } => left.is_mirror_of(right),
            }
        }

        let balanced_symbetric_trees: Vec<Tree<char>> = balanced_trees(5, 'x')
            .into_iter()
            .filter(|tree| is_symmetric(tree))
            .collect();
    }

    #[test]
    fn p59() {
        fn hbal_trees<T: fmt::Display + Copy>(height: usize, value: T) -> Vec<Tree<T>> {
            if height == 0 {
                vec![Tree::End]
            } else if height == 1 {
                vec![Tree::leaf(value)]
            } else {
                let subtrees_height_ng_1 = hbal_trees(height - 1, value);
                let subtrees_height_ng_2 = hbal_trees(height - 2, value);
                let mut res = vec![];
                res.extend(subtrees_height_ng_1.iter().flat_map(|tree_ng_1| {
                    subtrees_height_ng_2.iter().flat_map(move |tree_ng_2| {
                        vec![
                            Tree::node(value, tree_ng_1.clone(), tree_ng_2.clone()),
                            Tree::node(value, tree_ng_2.clone(), tree_ng_1.clone()),
                        ]
                    })
                }));
                res.extend(subtrees_height_ng_1.iter().flat_map(|tree_ng_1| {
                    subtrees_height_ng_1.iter().map(move |tree_ng_1_| {
                        Tree::node(value, tree_ng_1.clone(), tree_ng_1_.clone())
                    })
                }));
                res
            }
        }

        let trees = hbal_trees(3, 'x');
        assert_eq!(trees.len(), 15);
    }

    #[test]
    fn p60() {
        fn hbal_trees<T: fmt::Display + Copy>(height: usize, value: T) -> Vec<Tree<T>> {
            if height == 0 {
                vec![Tree::End]
            } else if height == 1 {
                vec![Tree::leaf(value)]
            } else {
                let subtrees_height_ng_1 = hbal_trees(height - 1, value);
                let subtrees_height_ng_2 = hbal_trees(height - 2, value);
                let mut res = vec![];
                res.extend(subtrees_height_ng_1.iter().flat_map(|tree_ng_1| {
                    subtrees_height_ng_2.iter().flat_map(move |tree_ng_2| {
                        vec![
                            Tree::node(value, tree_ng_1.clone(), tree_ng_2.clone()),
                            Tree::node(value, tree_ng_2.clone(), tree_ng_1.clone()),
                        ]
                    })
                }));
                res.extend(subtrees_height_ng_1.iter().flat_map(|tree_ng_1| {
                    subtrees_height_ng_1.iter().map(move |tree_ng_1_| {
                        Tree::node(value, tree_ng_1.clone(), tree_ng_1_.clone())
                    })
                }));
                res
            }
        }

        fn min_hbal_nodes(height: usize) -> usize {
            let trees = hbal_trees(3, 'x');
            trees
                .iter()
                .reduce(|t1, t2| {
                    if t1.node_count() < t2.node_count() {
                        t1
                    } else {
                        t2
                    }
                })
                .unwrap()
                .node_count()
        }

        assert_eq!(min_hbal_nodes(3), 4);
    }

    #[test]
    fn p61() {
        let tree = Tree::node('x', Tree::leaf('x'), Tree::end());
        assert_eq!(tree.leaf_count(), 1);
    }

    #[test]
    fn p62() {
        fn leaf_list<T: fmt::Display + Copy>(tree: &Tree<T>) -> Vec<T> {
            match tree {
                Tree::End => {
                    vec![]
                }
                Tree::Node { value, left, right } => match (left.as_ref(), right.as_ref()) {
                    (Tree::End, Tree::End) => {
                        vec![*value]
                    }
                    _ => {
                        let mut res = vec![];
                        res.extend(leaf_list(left));
                        res.extend(leaf_list(right));
                        res
                    }
                },
            }
        }
        let tree = Tree::node(
            'a',
            Tree::leaf('b'),
            Tree::node('c', Tree::leaf('d'), Tree::leaf('e')),
        );

        assert_eq!(leaf_list(&tree), vec!['b', 'd', 'e']);
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum Tree<T: fmt::Display> {
    Node {
        value: T,
        left: Box<Tree<T>>,
        right: Box<Tree<T>>,
    },
    End,
}

impl<T: fmt::Display> Tree<T> {
    fn node_count(&self) -> usize {
        match self {
            Tree::End => 0,
            Tree::Node {
                value: _,
                left,
                right,
            } => left.node_count() + right.node_count() + 1,
        }
    }
    fn leaf_count(&self) -> usize {
        match self {
            Tree::End => 0,
            Tree::Node {
                value: _,
                left,
                right,
            } => match (left.as_ref(), right.as_ref()) {
                (Tree::End, Tree::End) => 1,
                _ => left.leaf_count() + right.leaf_count(),
            },
        }
    }
    fn is_mirror_of(&self, other: &Tree<T>) -> bool {
        match &self {
            Tree::End => match &other {
                Tree::End => true,
                _ => false,
            },
            Tree::Node {
                value: _,
                left,
                right,
            } => match &other {
                Tree::End => false,
                Tree::Node {
                    value: _,
                    left: other_left,
                    right: other_right,
                } => left.is_mirror_of(other_right) && right.is_mirror_of(other_left),
            },
        }
    }
}

impl<T: fmt::Display> Tree<T> {
    fn node(value: T, left: Tree<T>, right: Tree<T>) -> Tree<T> {
        Tree::Node {
            value,
            left: Box::new(left),
            right: Box::new(right),
        }
    }

    fn leaf(value: T) -> Tree<T> {
        Tree::Node {
            value,
            left: Box::new(Tree::End),
            right: Box::new(Tree::End),
        }
    }

    fn end() -> Tree<T> {
        Tree::End
    }
}

impl<T: fmt::Display> fmt::Display for Tree<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match &self {
            Tree::Node { value, left, right } => write!(f, "T({} {} {})", value, left, right),
            Tree::End => write!(f, "."),
        }
    }
}
