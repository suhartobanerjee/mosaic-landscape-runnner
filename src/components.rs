use std::any::Any;

use bevy::input::keyboard::KeyboardInput;
use bevy::math::f32;
use bevy::scene::ron::de::Position;
use bevy::window::PrimaryWindow;
use bevy::prelude::*;
use bevy_rapier2d::na::iter;
use bevy_rapier2d::prelude::*;


pub const PLAYER_SPEED: f32 = 500.0;
pub const PLAYER_JUMP_IMPULSE: f32 = 20.0;
pub const BASE_ROAD: (f32, f32, f32)= (0.0, 70.0, 0.0);
pub const BASE_PLAYER: (f32, f32, f32)= (80.0, BASE_ROAD.1 + 60.0, 0.0);
pub const PLAYER_SIZE: f32 = 70.0;


#[derive(Component)]
pub struct Player {
    //position: Box<Landscape>,
    speed: f32,
    jump_impulse: f32,
    is_jumping: bool
}

pub fn spawn_player(mut commands: Commands,
                    asset_server: Res<AssetServer>) {

    commands.spawn(
        SpriteBundle{
            transform: Transform::from_xyz(BASE_PLAYER.0, BASE_PLAYER.1, BASE_PLAYER.2),
            texture: asset_server.load("sprites/character_maleAdventurer_run0.png"),
            ..Default::default()
        },
        )
        .insert(Player { 
            //position: Box::new(*landscape),
            speed: PLAYER_SPEED,
            jump_impulse: PLAYER_JUMP_IMPULSE,
            is_jumping: false
        })
        .insert(RigidBody::Dynamic)
        .insert(Collider::cuboid(0.1, 0.1))
        .insert(ColliderMassProperties::Mass(80.0))
        .insert(GravityScale(100.0))
        .insert(ActiveEvents::COLLISION_EVENTS)
        .insert(Velocity { linvel: Vec2::new(0.0, 0.0), angvel: 0.0 })
        .insert(LockedAxes::ROTATION_LOCKED);


}

#[derive(Component)]
pub struct GameCamera;

pub fn spawn_camera(mut commands: Commands,
                    window_query: Query<&Window, With<PrimaryWindow>>) {
    
    let window = window_query.get_single().unwrap();

    commands.spawn((
        Camera2dBundle {
            transform: Transform::from_xyz(window.width() / 2.0, window.height() / 2.0, BASE_PLAYER.2),
            ..Default::default()
        },
        GameCamera {},
        )
    );
}

pub fn player_movement(keyboard_input: Res<ButtonInput<KeyCode>>,
                       mut player_query: Query<(&mut Velocity, &mut Player), With<Player>>,
                       time: Res<Time>) {

    for (mut velocity, mut player) in player_query.iter_mut() {
        let mut direction = Vec2::ZERO;

        if keyboard_input.pressed(KeyCode::ArrowRight) {
            direction = Vec2::new(1.0, 0.0);
        }
        if keyboard_input.pressed(KeyCode::ArrowLeft) {
            direction = Vec2::new(-1.0, 0.0);
        }
        if keyboard_input.pressed(KeyCode::ArrowUp) && !player.is_jumping {
            direction = Vec2::new(0.0, player.jump_impulse);
            player.is_jumping = true;
        }

        velocity.linvel = direction * player.speed;
    }
}

pub fn update_bin_player(player_query: Query<(&Transform, &mut Player), With<Player>>) {
    
    for (position, player) in player_query.iter() {
        let current_bin = (position.translation.x / 120.0).ceil();

    }
}



pub fn set_jumping_false_if_touching_floor(mut contact_events: EventReader<CollisionEvent>,
                                           mut query: Query<(Entity, &mut Player)>) {
    for event in contact_events.read() {
        for (entity_player, mut player) in query.iter_mut() {
            if let CollisionEvent::Started(h1, h2, _) = event {
                if h1 == &entity_player || h2 == &entity_player {
                    player.is_jumping = false;
                }
            }
        }
    }

}

pub fn camera_movement(keyboard_input: Res<ButtonInput<KeyCode>>,
                       player_query: Query<&Velocity, With<Player>>,
                       mut camera_query: Query<&mut Transform, With<GameCamera>>,
                       time: Res<Time>) {

    if let Ok(mut transform) = camera_query.get_single_mut() {
        let mut direction = Vec3::ZERO;
        let player_x = player_query.get_single().unwrap().linvel.x;

        if keyboard_input.pressed(KeyCode::ArrowRight) {
            direction = Vec3::new(1.0, 0.0, 0.0)
        }
        if keyboard_input.pressed(KeyCode::ArrowLeft) {
            direction = Vec3::new(-1.0, 0.0, 0.0)
        }

        if direction.length() > 0.0 {
            direction = direction.normalize();
        }

        transform.translation += Vec3::new(player_x, 0.0, 0.0) * time.delta_seconds() * 0.5;
    }
}



pub fn confine_player_movement(mut player_query: Query<&mut Transform, With<Player>>,
                               window_query: Query<&Window, With<PrimaryWindow>>) {

    if let Ok(mut player_transform) = player_query.get_single_mut() {
        
        let window = window_query.get_single().unwrap();
        let x_min = BASE_PLAYER.0;
        let x_max = window.width() - PLAYER_SIZE / 2.0;
        let y_min = BASE_PLAYER.1;
        let y_max = y_min + BASE_ROAD.1;

        let mut translation = player_transform.translation;

        //bound the movement
        if translation.x < x_min {
           translation.x = x_min; 
        }
        if translation.x > x_max {
           translation.x = x_max; 
        }
        if translation.y > y_max {
           translation.y = y_max; 
        }
        if translation.y < y_min {
           translation.y = y_min; 
        }

        player_transform.translation = translation;
    }
    
}


pub fn confine_camera_movement(mut camera_query: Query<&mut Transform, With<GameCamera>>,
                               window_query: Query<&Window, With<PrimaryWindow>>) {

    if let Ok(mut camera_transform) = camera_query.get_single_mut() {
        
        let window = window_query.get_single().unwrap();
        let x_min = BASE_PLAYER.0;
        let x_max = window.width() - PLAYER_SIZE / 2.0;
        let y_min = BASE_PLAYER.1;
        let y_max = y_min + BASE_ROAD.1;

        let mut translation = camera_transform.translation;

        //bound the movement
        if translation.x < x_min {
           translation.x = x_min; 
        }
        if translation.x > x_max {
           translation.x = x_max; 
        }
        if translation.y > y_max {
           translation.y = y_max; 
        }
        if translation.y < y_min {
           translation.y = y_min; 
        }

        camera_transform.translation = translation;
    }
    
}

#[derive(Debug, Clone, Copy)]
pub enum SvState {
    Ref,
    Amp,
    Del,
    Inv
}


#[derive(Copy, Clone, Debug)]
pub struct Land {
    bin_id: u32,
    sv_state: SvState,
}

#[derive(Resource, Clone, Debug)]
pub struct Landscape {
    map: Vec<Land>
}

impl Landscape {
    fn new() -> Self {
        return Self { map: Vec::new() };
    }
}


pub fn spawn_landscape(mut commands: Commands,
                       asset_server: Res<AssetServer>) {
    

    let offset = 120.0;
    let rounds: u16 = 10;
    let mut sprite_dir: &str;
    let mut bins: u32;
    let mut sv_type: SvState;
    let mut landscape: Landscape = Landscape::new();


    for round in 0..rounds{

        if round < 5 {
            sv_type = SvState::Ref;
            sprite_dir = "sprites/roadTexture_13.png";
        } else {
            sv_type = SvState::Del;
            sprite_dir = "sprites/roadTexture_85.png";
        }


        let sprite_bundle = SpriteBundle{
            transform: Transform::from_xyz(BASE_ROAD.0 + offset * f32::from(round), BASE_ROAD.1, BASE_ROAD.2 - 1.0),
            texture: asset_server.load(sprite_dir),
            ..Default::default()
        };

        landscape.map.push(Land { bin_id: round as u32, sv_state: sv_type });

        commands.insert_resource(landscape.clone());

        commands.spawn(sprite_bundle)
        .insert(RigidBody::Fixed)
        .insert(Collider::cuboid(70.0, 70.0))
        .insert(ActiveEvents::COLLISION_EVENTS);
        
    }
}

pub fn check_sv(player_query: Query<(&Transform, &Player)>,
                landscape_query: Res<Landscape>) {


    let player_position = player_query.get_single().expect("player_position not found").0.translation;
    let current_bin = (player_position.x / 100.0).floor() as usize;
    let sv_state = landscape_query.map[current_bin];

    //println!("player: {}", player_position.x);
    // println!("player: {:?}", player_query.get_single().unwrap().1.position);
    println!("land: {:?}", sv_state);
    // println!("land_type: {:?}", landscape_instance.1.sv_state);
}


